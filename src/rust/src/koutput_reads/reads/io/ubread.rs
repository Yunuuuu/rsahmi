use std::fs::File;
use std::io::BufWriter;
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Error, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use indicatif::{MultiProgress, ProgressBar, ProgressFinish};
use rustc_hash::FxHashSet as HashSet;

use crate::batchsender::BatchSender;
use crate::koutput_reads::range::*;
use crate::kractor::reads::io::reader::BytesChunkPairedReader;
use crate::kractor::reads::parser::fasta::FastaRecordWithUMIBarcode;
use crate::kractor::reads::parser::fastq::FastqContainer;

pub fn reader_kractor_ubread_read(
    id_sets: HashSet<&[u8]>,
    fq: &str,
    ofile: &str,
    ubread: &str,
    umi_ranges: Vec<RangeKind>,
    barcode_ranges: Vec<RangeKind>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Create output file and wrap in buffered writer
    let mut writer = BufWriter::with_capacity(buffer_size, File::create(ofile)?);

    // Open and memory-map the input FASTQ file
    let reader1 = File::open(fq)?;
    let reader2 = File::open(ubread)?;
    let size1 = reader1.metadata()?.len();
    let size2 = reader2.metadata()?.len();
    let style = crate::progress_style()?;

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<FastaRecordWithUMIBarcode<Bytes>>>,
            Receiver<Vec<FastaRecordWithUMIBarcode<Bytes>>>,
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
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let m = MultiProgress::new();
        let pb1 = m.add(ProgressBar::new(size1).with_finish(ProgressFinish::Abandon));
        pb1.set_prefix("Parsing reads");
        pb1.set_style(style.clone());
        let pb2 = m.add(ProgressBar::new(size2).with_finish(ProgressFinish::Abandon));
        pb2.set_prefix("Parsing ubread");
        pb2.set_style(style);
        let mut reader = BytesChunkPairedReader::with_capacity(chunk_size, reader1, reader2);
        reader.set_label1("reads");
        reader.set_label2("ubread");

        #[cfg(not(test))]
        reader.attach_bars(pb1, pb2);

        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `id_sets`
            // Shared atomic flag to signal if any thread encountered an error
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx): (Sender<Error>, Receiver<Error>) = crossbeam_channel::bounded(1);
            // Reuse the shared ID set for matching records
            let id_sets = &id_sets;

            // Reuse the umi_ranges and barcode_ranges
            let umi_ranges = &umi_ranges;
            let barcode_ranges = &barcode_ranges;
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
                    let thread_err_tx = err_tx.clone();
                    s.spawn(move |_| {
                        let mut container1 = FastqContainer::new();
                        let mut container2 = FastqContainer::new();
                        loop {
                            match chunk.read_record(&mut container1, &mut container2) {
                                Ok(value) => match value {
                                    Some((record1, record2)) => {
                                        if id_sets.contains(&record1.id[..]) {
                                            let umi_result = extract_pattern_from_sequence(
                                                &record2.seq[..],
                                                umi_ranges,
                                            );
                                            let umi = match umi_result {
                                                Ok(umi) => umi,
                                                Err(e) => {
                                                    thread_has_error.store(true, Relaxed);
                                                    let _ = thread_err_tx.try_send(e);
                                                    return;
                                                }
                                            };
                                            let barcode_result = extract_pattern_from_sequence(
                                                &record2.seq[..],
                                                barcode_ranges,
                                            );
                                            // Try extracting barcode
                                            let barcode = match barcode_result {
                                                Ok(barcode) => barcode,
                                                Err(e) => {
                                                    thread_has_error.store(true, Relaxed);
                                                    let _ = thread_err_tx.try_send(e);
                                                    return;
                                                }
                                            };
                                            let record = FastaRecordWithUMIBarcode::new(
                                                record1, umi, barcode,
                                            );
                                            // Attempt to send the matching record to the writer thread.
                                            match thread_tx.send(record) {
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
