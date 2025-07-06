use std::fs::File;
use std::io::{BufWriter, Read};
use std::iter::zip;
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Error, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::batchsender::BatchSender;
use crate::parser::fastq::FastqRecord;
use crate::{fastq_reader, seq_action::*};

pub(crate) fn reader_seq_refine_paired_read<R1: Read + Send, R2: Read + Send>(
    reader1: R1,
    ofile1: Option<&str>,
    reader2: R2,
    ofile2: Option<&str>,
    ref actions: SubseqPairedActions,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Create output file and wrap in buffered writer
    let writer1;
    if let Some(file) = ofile1 {
        let file = File::create(&file)?;
        let bw = BufWriter::with_capacity(buffer_size, file);
        writer1 = Some(GzEncoder::new(bw, Compression::new(4)));
    } else {
        writer1 = None
    }
    let writer2;
    if let Some(file) = ofile2 {
        let file = File::create(&file)?;
        let bw = BufWriter::with_capacity(buffer_size, file);
        writer2 = Some(GzEncoder::new(bw, Compression::new(4)));
    } else {
        writer2 = None
    }

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<(FastqRecord<Bytes>, FastqRecord<Bytes>)>>,
            Receiver<Vec<(FastqRecord<Bytes>, FastqRecord<Bytes>)>>,
        ) = crate::new_channel(nqueue);

        let (writer1_tx, writer1_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let (writer2_tx, writer2_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let (reader1_tx, reader1_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let (reader2_tx, reader2_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let has_writer1 = writer1.is_some();
        let writer1_handle = if let Some(mut writer) = writer1 {
            Some(scope.spawn(move || -> Result<()> {
                for chunk in writer1_rx {
                    for record in chunk {
                        record.write(&mut writer)?;
                    }
                }
                Ok(())
            }))
        } else {
            None
        };

        let has_writer2 = writer2.is_some();
        let writer2_handle = if let Some(mut writer) = writer2 {
            Some(scope.spawn(move || -> Result<()> {
                for chunk in writer2_rx {
                    for record in chunk {
                        record.write(&mut writer)?;
                    }
                }
                Ok(())
            }))
        } else {
            None
        };

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            let mut writer1_tx = BatchSender::with_capacity(batch_size, writer1_tx);
            let mut writer2_tx = BatchSender::with_capacity(batch_size, writer2_tx);
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for (record1, record2) in chunk {
                    if has_writer1 {
                        writer1_tx.send(record1)?;
                    }
                    if has_writer2 {
                        writer2_tx.send(record2)?;
                    }
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `id_sets`
            // Shared atomic flag to signal if any thread encountered an error
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            // Rayon scope for spawning parsing threads
            rayon::scope(|s| -> Result<()> {
                loop {
                    let (records1, records2) = match (reader1_rx.recv(), reader2_rx.recv()) {
                        (Ok(rec1), Ok(rec2)) => (rec1, rec2),
                        (Err(_), Ok(_)) => {
                            return Err(anyhow!(
                                "FASTQ pairing error: read1 channel closed unexpectedly before read2"
                            ));
                        }
                        (Ok(_), Err(_)) => {
                            return Err(anyhow!(
                                "FASTQ pairing error: read2 channel closed unexpectedly before read1"
                            ));
                        }
                        (Err(_), Err(_)) => {
                            return Ok(());
                        }
                    };

                    if records1.len() != records2.len() {
                        return Err(anyhow!(
                            "FASTQ pairing error: mismatched record counts (read1: {}, read2: {})",
                            records1.len(),
                            records2.len()
                        ));
                    }

                    // If an error already occurred, stop spawning new threads
                    if has_error.load(Relaxed) {
                        return Ok(());
                    }

                    // Initialize a thread-local batch sender for matching records
                    let mut thread_tx = BatchSender::with_capacity(batch_size, parser_tx.clone());

                    // Clone the shared error state for this thread
                    let thread_has_error = has_error.clone();
                    let thread_err_tx: Sender<Error> = err_tx.clone();
                    s.spawn(move |_| {
                        for (mut record1, mut record2) in zip(records1, records2) {
                            match actions.transform_fastq_bytes(&mut record1, &mut record2) {
                                Ok(_) => {
                                    // Attempt to send the matching record to the writer thread.
                                    match thread_tx.send((record1, record2)) {
                                        Ok(_) => continue,
                                        Err(_) => {
                                            // If this fails, it means the writer thread has already exited due to an error.
                                            // Since that error will be reported separately, we can safely ignore this send failure.
                                            // Mark that an error has occurred
                                            thread_has_error.store(true, Relaxed);
                                            return ();
                                        }
                                    };
                                }
                                Err(e) => {
                                    // Mark that an error has occurred
                                    thread_has_error.store(true, Relaxed);
                                    // Only the first error is reported
                                    let _ = thread_err_tx.try_send(e.into());
                                    return ();
                                }
                            }
                        }
                        let _ = thread_tx.flush();
                    });
                }
            })?;
            // Clean up: close the error channel and drop parser sender
            drop(err_tx);
            drop(parser_tx);

            // Report the first error if any thread encountered one
            if has_error.load(Relaxed) {
                let err = err_rx.recv()?; // Safe unwrap because has_error is true
                Err(anyhow!(err))
            } else {
                Ok(())
            }
        });

        // ─── reader Thread ─────────────────────────────────────
        let reader1_handle = scope.spawn(move || -> Result<()> {
            let mut reader1 = fastq_reader::FastqReader::new(reader1);
            let mut reader1_tx = BatchSender::with_capacity(chunk_size, reader1_tx);
            while let Some(record) = reader1.read_record()? {
                reader1_tx.send(record)?;
            }
            reader1_tx.flush()?;
            Ok(())
        });

        let reader2_handle = scope.spawn(move || -> Result<()> {
            let mut reader2 = fastq_reader::FastqReader::new(reader2);
            let mut reader2_tx = BatchSender::with_capacity(chunk_size, reader2_tx);
            while let Some(record) = reader2.read_record()? {
                reader2_tx.send(record)?;
            }
            reader2_tx.flush()?;
            Ok(())
        });

        // ─── Join Threads and Propagate Errors ────────────────
        if let Some(writer_handle) = writer1_handle {
            writer_handle
                .join()
                .map_err(|e| anyhow!("Writer1 thread panicked: {:?}", e))??
        };
        if let Some(writer_handle) = writer2_handle {
            writer_handle
                .join()
                .map_err(|e| anyhow!("Writer2 thread panicked: {:?}", e))??
        };
        writer_handle
            .join()
            .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))??;
        parser_handle
            .join()
            .map_err(|e| anyhow!("Parser thread panicked: {:?}", e))??;
        reader1_handle
            .join()
            .map_err(|e| anyhow!("Writer2 thread panicked: {:?}", e))??;
        reader2_handle
            .join()
            .map_err(|e| anyhow!("Writer2 thread panicked: {:?}", e))??;
        Ok(())
    })
}
