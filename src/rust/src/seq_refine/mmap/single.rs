use std::fs::File;
use std::io::BufWriter;
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Error, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};

use crate::batchsender::BatchSender;
use crate::parser::fastq::FastqRecord;
use crate::parser::fastq::FastqSliceChunkReader;
use crate::reader::slice::SliceProgressBarReader;
use crate::seq_action::*;

pub fn mmap_seq_refine_single_read(
    reader: SliceProgressBarReader,
    ofile: &str,
    ref actions: SubseqActions, // ✅ Use your unified struct
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Open output file and wrap in buffered writer
    let mut writer = BufWriter::with_capacity(buffer_size, File::create(ofile)?);

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<(FastqRecord<&[u8]>, (Vec<Bytes>, Bytes, Bytes))>>,
            Receiver<Vec<(FastqRecord<&[u8]>, (Vec<Bytes>, Bytes, Bytes))>>,
        ) = crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for (mut record, (tags, seq, qual)) in chunk {
                    let desc = record.desc.map_or_else(
                        || tags.iter().flatten().copied().collect::<Bytes>(),
                        |d| {
                            d.as_ref()
                                .into_iter()
                                .chain(std::iter::once(&b' '))
                                .chain(tags.iter().flatten())
                                .copied()
                                .collect::<Bytes>()
                        },
                    );
                    record.desc = Some(&desc);
                    record.seq = &seq;
                    record.qual = &qual;
                    record.write(&mut writer)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let mut reader = FastqSliceChunkReader::with_capacity(chunk_size, reader);

        let parser_handle = scope.spawn(move || -> Result<()> {
            // will move `reader`, `parser_tx`, and `id_sets`
            // Shared atomic flag to signal if any thread encountered an error
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            // Reuse the shared ID set for matching records
            // Rayon scope for spawning parsing threads
            rayon::scope(|s| {
                while let Some(mut chunk) = reader.chunk_reader() {
                    // If an error already occurred, stop spawning new threads
                    if has_error.load(Relaxed) {
                        return ();
                    }
                    // Initialize a thread-local batch sender for matching records
                    let mut thread_tx = BatchSender::with_capacity(batch_size, parser_tx.clone());
                    // Clone the shared error state for this thread
                    let thread_has_error = has_error.clone();
                    let thread_err_tx: Sender<Error> = err_tx.clone();
                    s.spawn(move |_| {
                        loop {
                            match chunk.read_record() {
                                Ok(opt_record) => match opt_record {
                                    Some(record) => {
                                        match actions.apply_to_fastq(&record) {
                                            Ok((tags, seq, qual)) => {
                                                // Attempt to send the matching record to the writer thread.
                                                match thread_tx.send((record, (tags, seq, qual))) {
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
            });
            // Clean up: close the error channel and drop parser sender
            drop(err_tx);
            drop(parser_tx);

            // Report the first error if any thread encountered one
            if has_error.load(Relaxed) {
                let err = err_rx.recv()?; // Safe unwrap because has_error is true
                Err(err)
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
