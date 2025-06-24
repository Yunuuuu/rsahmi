use std::fs::File;
use std::io::BufWriter;
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Result};
use crossbeam_channel::{Receiver, Sender};
#[cfg(unix)]
use memmap2::Advice;
use memmap2::Mmap;
use rustc_hash::FxHashSet as HashSet;

use super::reader::SliceChunkPairedReader;
use crate::batchsender::BatchSender;
use crate::parser::fasta::FastaRecord;

pub fn mmap_kractor_paired_read(
    id_sets: HashSet<&[u8]>,
    fq1: &str,
    ofile1: &str,
    fq2: &str,
    ofile2: &str,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Create output file and wrap in buffered writer
    let mut writer1 =
        BufWriter::with_capacity(buffer_size, File::create(ofile1)?);
    let mut writer2 =
        BufWriter::with_capacity(buffer_size, File::create(ofile2)?);

    // Open and memory-map the input FASTQ file
    let file1 = File::open(fq1)?;
    let map1 = unsafe { Mmap::map(&file1) }?;
    // #[cfg(unix)]
    // map1.advise(Advice::Sequential)?;

    let file2 = File::open(fq2)?;
    let map2 = unsafe { Mmap::map(&file2) }?;
    // #[cfg(unix)]
    // map2.advise(Advice::Sequential)?;

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<(FastaRecord<&[u8]>, FastaRecord<&[u8]>)>>,
            Receiver<Vec<(FastaRecord<&[u8]>, FastaRecord<&[u8]>)>>,
        ) = crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for (record1, record2) in chunk {
                    record1.write(&mut writer1)?;
                    record2.write(&mut writer2)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let mut reader =
            SliceChunkPairedReader::with_capacity(chunk_size, &map1, &map2);
        reader.set_label1("fq1");
        reader.set_label2("fq2");

        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `matcher`
            // Shared atomic flag to signal if any thread encountered an error
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            // Reuse the shared ID set for matching records
            let id_sets = &id_sets;
            // Rayon scope for spawning parsing threads
            rayon::scope(|s| -> Result<()> {
                for parser_result in reader {
                    let mut parser = parser_result?;
                    // If an error already occurred, stop spawning new threads
                    if has_error.load(Relaxed) {
                        return Ok(());
                    }
                    // Initialize a thread-local batch sender for matching records
                    let mut thread_tx = BatchSender::with_capacity(
                        batch_size,
                        parser_tx.clone(),
                    );
                    // Clone the shared error state for this thread
                    let thread_has_error = has_error.clone();
                    let thread_err_tx = err_tx.clone();
                    s.spawn(move |_| {
                        loop {
                            match parser.read_record() {
                                Ok(value) => match value {
                                    Some((record1, record2)) => {
                                        if id_sets.contains(record1.id) {
                                            // Attempt to send the matching record to the writer thread.
                                            // If this fails, it means the writer thread has already exited due to an error.
                                            // Since that error will be reported separately, we can safely ignore this send failure.
                                            match thread_tx
                                                .send((record1, record2))
                                            {
                                                Ok(_) => continue,
                                                Err(_) => return (),
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
