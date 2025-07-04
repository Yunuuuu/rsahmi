use std::fs::File;
use std::io::{BufWriter, Read};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use rustc_hash::FxHashSet as HashSet;

use crate::batchsender::BatchSender;
use crate::parser::fasta::FastaRecord;
use crate::parser::fastq::FastqBytesChunkReader;
use crate::parser::fastq::FastqContainer;
use crate::reader::bytes::BytesProgressBarReader;

/// Reads a single-end FASTQ file in chunks, filters records by ID set, and writes matching records to output.
/// This function uses thread + rayon + channel architecture for parallelism.
pub fn reader_kractor_single_read<R: Read + Send>(
    id_sets: HashSet<&[u8]>,
    reader: BytesProgressBarReader<R>,
    ofile: &str,
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
            Sender<Vec<FastaRecord<Bytes>>>,
            Receiver<Vec<FastaRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for record in chunk {
                    record.write(&mut writer)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Reads and splits the input FASTQ into chunks.
        // Each chunk is parsed in parallel using Rayon.
        // Matching records (based on `id_sets`) are sent to the writer via the channel.
        let mut reader = FastqBytesChunkReader::with_capacity(chunk_size, reader);

        let parser_handle = scope.spawn(move || -> Result<()> {
            // will move `reader`, `parser_tx`, and `id_sets`
            // Flag if any thread fails
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            // Reuse the shared ID set for matching records
            let id_sets = &id_sets;
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
                    // Only restore the parsing error, for sending error, it means writer thread has already exited due to an error.
                    let thread_err_tx = err_tx.clone();
                    s.spawn(move |_| {
                        let mut container = FastqContainer::new();
                        loop {
                            match chunk.read_record(&mut container) {
                                Ok(value) => match value {
                                    Some(record) => {
                                        if id_sets.contains(&record.id[..]) {
                                            // Attempt to send the matching record to the writer thread.
                                            match thread_tx.send(record.into()) {
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
                                    }
                                    None => break,
                                },
                                Err(e) => {
                                    // Mark that an error has occurred
                                    thread_has_error.store(true, Relaxed);
                                    // Only the first error is reported
                                    let _ = thread_err_tx.try_send(e);
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
