use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::{CompressionLvl, Compressor};
use rustc_hash::FxHashSet as HashSet;

use crate::batchsender::BatchSender;
use crate::fastq_reader::*;
use crate::fastq_record::FastqRecord;
use crate::utils::*;

pub(super) fn parse_single<P: AsRef<Path> + ?Sized>(
    id_sets: &HashSet<&[u8]>,
    input_path: &P,
    input_bar: Option<ProgressBar>,
    output_path: &P,
    output_bar: Option<ProgressBar>,
    compression_level: i32,
    batch_size: usize,
    chunk_bytes: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let input: &Path = input_path.as_ref();
    let output: &Path = output_path.as_ref();

    // Ensure compression level is validated and converted before entering thread scope.
    // Doing this outside avoids redundant validation across parser threads.
    let compression_level = CompressionLvl::new(compression_level)
        .map_err(|e| anyhow!("Invalid 'compression_level': {:?}", e))?;
    std::thread::scope(|scope| -> Result<()> {
        // Two communication pipelines are set up to decouple IO and CPU-intensive work:
        // - reader_tx: transfers raw FASTQ records to parser threads
        // - writer_tx: receives compressed byte chunks from parser threads
        let (writer_tx, writer_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = new_channel(nqueue);
        let (reader_tx, reader_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // A single thread handles file output to ensure atomic write order and leverage buffered IO.
        // This thread consumes compressed chunks, not raw records, for performance.
        let writer_handle = scope.spawn(move || -> Result<()> {
            let mut writer = BufWriter::with_capacity(chunk_bytes, new_writer(output, output_bar)?);

            // Iterate over each received batch of records
            for chunk in writer_rx {
                writer
                    .write_all(&chunk)
                    .with_context(|| format!("(Writer) Failed to write FastqRecord to output"))?;
            }
            writer
                .flush()
                .with_context(|| format!("(Writer) Failed to flush writer"))?;
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Multiple parser threads are used to exploit CPU parallelism for `transform_fastq` and gzip compression.
        // Each thread transforms records and buffers them into a local pool,
        // which is periodically flushed into the writer pipeline.
        let mut parser_handles = Vec::with_capacity(threads);
        let gzip = gz_compressed(output);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                // Temporary buffer for current output chunk
                let mut records_pool: Vec<u8> = Vec::with_capacity(chunk_bytes);
                let mut compressor = Compressor::new(compression_level);
                while let Ok(records) = rx.recv() {
                    for record in records {
                        if id_sets.contains(record.id.as_ref()) {
                            // Flush when pool is too full to accept the next record.
                            // This ensures output chunks remain near the target block size.
                            if records_pool.capacity() - records_pool.len() < record.bytes_size() {
                                let mut pack = Vec::with_capacity(chunk_bytes);
                                std::mem::swap(&mut records_pool, &mut pack);
                                // Compress if gzip file
                                if gzip {
                                    pack = gzip_pack(&pack, &mut compressor)?
                                }

                                // Send compressed or raw bytes to writer
                                tx.send(pack).with_context(|| {
                                    format!(
                                        "(Parser) Failed to send parsed record to Writer thread"
                                    )
                                })?;
                            }
                        }
                        // Append encoded record to buffer
                        record.extend(&mut records_pool);
                    }
                }

                // Flush remaining records if any
                if !records_pool.is_empty() {
                    let pack = if gzip {
                        gzip_pack(&records_pool, &mut compressor)?
                    } else {
                        records_pool
                    };
                    tx.send(pack).with_context(|| {
                        format!("(Parser) Failed to send parsed record to Writer thread")
                    })?;
                }
                Ok(())
            });
            parser_handles.push(handle);
        }
        drop(reader_rx);
        drop(writer_tx);

        // ─── reader Thread ─────────────────────────────────────
        let reader_handle = scope.spawn(move || -> Result<()> {
            let mut reader =
                FastqReader::with_capacity(BUFFER_SIZE, new_reader(input, BUFFER_SIZE, input_bar)?);
            let mut reader_tx = BatchSender::with_capacity(batch_size, reader_tx);
            while let Some(record) = reader
                .read_record()
                .with_context(|| format!("(Reader) Failed to read FASTQ record"))?
            {
                reader_tx.send(record).with_context(|| {
                    format!("(Reader) Failed to send FASTQ records to Parser thread")
                })?;
            }
            reader_tx.flush().with_context(|| {
                format!("(Reader) Failed to flush FASTQ records to Parser thread")
            })?;
            Ok(())
        });

        // ─── Join Threads and Propagate Errors ────────────────
        writer_handle
            .join()
            .map_err(|e| anyhow!("(Writer) thread panicked: {:?}", e))??;
        for handler in parser_handles {
            handler
                .join()
                .map_err(|e| anyhow!("(Parser) thread panicked: {:?}", e))??;
        }
        reader_handle
            .join()
            .map_err(|e| anyhow!("(Reader) thread panicked: {:?}", e))??;
        Ok(())
    })
}
