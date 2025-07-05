use std::fs::File;
use std::io::{BufWriter, Read};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Error, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};

use crate::batchsender::BatchSender;
use crate::parser::fastq::{FastqBytesChunkPairedReader, FastqContainer, FastqRecord};
use crate::reader::bytes::BytesProgressBarReader;
use crate::seq_action::*;

pub(super) fn reader_seq_refine_paired_read<R: Read + Send>(
    reader1: BytesProgressBarReader<R>,
    ofile1: Option<&str>,
    reader2: BytesProgressBarReader<R>,
    ofile2: Option<&str>,
    ref actions: SubseqPairedActions,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Create output file and wrap in buffered writer
    let mut writer1;
    if let Some(file) = ofile1 {
        writer1 = Some(BufWriter::with_capacity(buffer_size, File::create(file)?));
    } else {
        writer1 = None
    }
    let mut writer2;
    if let Some(file) = ofile2 {
        writer2 = Some(BufWriter::with_capacity(buffer_size, File::create(file)?));
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

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for (record1, record2) in chunk {
                    if let Some(ref mut writer) = &mut writer1 {
                        record1.write(writer)?;
                    }
                    if let Some(ref mut writer) = &mut writer2 {
                        record2.write(writer)?;
                    }
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let mut reader = FastqBytesChunkPairedReader::with_capacity(chunk_size, reader1, reader2);

        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `id_sets`
            // Shared atomic flag to signal if any thread encountered an error
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            // Rayon scope for spawning parsing threads
            rayon::scope(|s| -> Result<()> {
                while let Some(chunk) = reader.chunk_reader()? {
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
                        let mut container1 = FastqContainer::new();
                        let mut container2 = FastqContainer::new();
                        loop {
                            match chunk.read_record(&mut container1, &mut container2) {
                                Ok(value) => match value {
                                    Some((mut record1, mut record2)) => {
                                        match actions
                                            .transform_fastq_bytes(&mut record1, &mut record2)
                                        {
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
                                    None => break,
                                },
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
                Ok(())
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

        // ─── Join Threads and Propagate Errors ────────────────
        writer_handle
            .join()
            .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))??;
        parser_handle
            .join()
            .map_err(|e| anyhow!("Parser thread panicked: {:?}", e))??;
        Ok(())
    })
}
