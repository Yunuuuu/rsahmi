use std::fs::File;
use std::io::BufWriter;

use anyhow::{anyhow, Result};
use crossbeam_channel::{Receiver, Sender};
use memmap2::{Advice, Mmap};
use rayon::prelude::*;
use rustc_hash::FxHashSet as HashSet;

use super::reader::SliceChunkPairedReader;
use crate::batchsender::BatchSender;
use crate::kractor::reads::range::*;
use crate::parser::fasta::FastaRecordWithUMIBarcode;

pub fn mmap_kractor_ubread_read(
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
    let mut writer =
        BufWriter::with_capacity(buffer_size, File::create(ofile)?);

    // Open and memory-map the input FASTQ file
    let file1 = File::open(fq)?;
    let map1 = unsafe { Mmap::map(&file1) }?;
    map1.advise(Advice::Sequential)?;

    let file2 = File::open(ubread)?;
    let map2 = unsafe { Mmap::map(&file2) }?;
    map2.advise(Advice::Sequential)?;

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<FastaRecordWithUMIBarcode<&[u8]>>>,
            Receiver<Vec<FastaRecordWithUMIBarcode<&[u8]>>>,
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
        let mut reader =
            SliceChunkPairedReader::with_capacity(chunk_size, &map1, &map2);
        reader.set_label1("reads");
        reader.set_label2("ubread");

        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `id_sets`
            reader.par_bridge().try_for_each_init(
                // Initialize per-thread batching sender
                || BatchSender::with_capacity(batch_size, parser_tx.clone()),
                |thread_tx, parser_result| -> Result<()> {
                    let mut parser = parser_result?;
                    // Each thread parses a chunk, filters it, and sends matching records
                    while let Some((record1, record2)) = parser.read_record()? {
                        if id_sets.contains(record1.id) {
                            let umi = extract_pattern_from_sequence(
                                record2.seq,
                                &umi_ranges,
                            )?;
                            let barcode = extract_pattern_from_sequence(
                                record2.seq,
                                &barcode_ranges,
                            )?;
                            let record = FastaRecordWithUMIBarcode::new(
                                record1, umi, barcode,
                            );
                            // Wrap `SendError` in `anyhow::Error` explicitly.
                            // this avoids capturing internal references (from `mmap`).
                            // whose lifetime is not long enough to hold it
                            thread_tx.send(record).map_err(|e| {
                                anyhow!(
                                    "Failed to send to Writer thread: {}",
                                    e
                                )
                            })?;
                        }
                    }
                    thread_tx.flush().map_err(|e| {
                        anyhow!("Failed to send to Writer thread: {}", e)
                    })?;
                    Ok(())
                },
            )
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
