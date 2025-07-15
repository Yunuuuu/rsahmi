use std::io::{BufWriter, Write};
use std::path::Path;

use aho_corasick::AhoCorasick;
use anyhow::{anyhow, Context, Result};
use bytes::{BufMut, BytesMut};
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::{CompressionLvl, Compressor};
use memchr::memchr;
use rustc_hash::FxHashSet as HashSet;

use crate::batchsender::BatchSender;
use crate::reader0::LineReader;
use crate::utils::*;

pub(super) fn parse_koutput<P: AsRef<Path> + ?Sized>(
    input_path: &P,
    input_bar: Option<ProgressBar>,
    output_path: &P,
    output_bar: Option<ProgressBar>,
    include_sets: HashSet<&[u8]>,
    exclude_aho: Option<AhoCorasick>,
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
        let (reader_tx, reader_rx): (Sender<Vec<BytesMut>>, Receiver<Vec<BytesMut>>) =
            new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // A single thread handles file output to ensure atomic write order and leverage buffered IO.
        // This thread consumes compressed chunks, not raw records, for performance.
        let writer_handle = scope.spawn(move || -> Result<()> {
            let mut writer = BufWriter::with_capacity(chunk_bytes, new_writer(output, output_bar)?);

            // Iterate over each received batch of records
            for chunk in writer_rx {
                writer
                    .write_all(&chunk)
                    .with_context(|| format!("(Writer) Failed to write Fastq records to output"))?;
            }
            writer
                .flush()
                .with_context(|| format!("(Writer) Failed to flush writer"))?;
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams Kraken2 output data, filters by ID set
        let mut parser_handles = Vec::with_capacity(threads);
        let gzip = gz_compressed(output);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let include_sets = &include_sets;
            let exclude_aho = &exclude_aho;
            let handle = scope.spawn(move || -> Result<()> {
                let mut pool: Vec<u8> = Vec::with_capacity(chunk_bytes);
                let mut compressor = Compressor::new(compression_level);
                while let Ok(lines) = rx.recv() {
                    for line in lines {
                        if kractor_match_aho(&include_sets, &exclude_aho, &line) {
                            // Flush when pool is too full to accept the next record.
                            // This ensures output chunks remain near the target block size.
                            if pool.capacity() - pool.len() < (line.len() + 1) {
                                let mut pack = Vec::with_capacity(chunk_bytes);
                                std::mem::swap(&mut pool, &mut pack);
                                // Compress if gzip file
                                if gzip {
                                    pack = gzip_pack(&pack, &mut compressor)?
                                }

                                // Send compressed or raw bytes to writer
                                tx.send(pack).with_context(|| {
                                    format!("(Parser) Failed to send parsed lines to Writer thread")
                                })?;
                            }
                            // Append encoded lines to buffer
                            pool.extend_from_slice(&line);
                            pool.put_u8(b'\n');
                        };
                    }
                }
                // Flush remaining lines if any
                if !pool.is_empty() {
                    let pack = if gzip {
                        gzip_pack(&pool, &mut compressor)?
                    } else {
                        pool
                    };
                    tx.send(pack).with_context(|| {
                        format!("(Parser) Failed to send parsed lines to Writer thread")
                    })?;
                };
                Ok(())
            });
            parser_handles.push(handle);
        }
        drop(reader_rx);
        drop(writer_tx);

        // ─── reader Thread ─────────────────────────────────────
        let reader_handle = scope.spawn(move || -> Result<()> {
            let mut reader =
                LineReader::with_capacity(BUFFER_SIZE, new_reader(input, BUFFER_SIZE, input_bar)?);
            let mut reader_tx = BatchSender::with_capacity(batch_size, reader_tx);
            while let Some(record) = reader
                .read_line()
                .with_context(|| format!("(Reader) Failed to read line"))?
            {
                reader_tx
                    .send(record)
                    .with_context(|| format!("(Reader) Failed to send lines to Parser thread"))?;
            }
            reader_tx
                .flush()
                .with_context(|| format!("(Reader) Failed to flush lines to Parser thread"))?;
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

fn kractor_match_aho(
    include_sets: &HashSet<&[u8]>,
    exclude_aho: &Option<AhoCorasick>,
    line: &[u8],
) -> bool {
    let mut field_start = 0usize;
    let mut field_index = 0usize;
    while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
        if field_index == 2 {
            let field = &line[field_start .. (field_start + tab_pos)];
            if let Some(start) = KOUTPUT_TAXID_PREFIX_FINDER.find(field) {
                let start = start + KOUTPUT_TAXID_PREFIX.len();
                if let Some(end) = memchr(KOUTPUT_TAXID_SUFFIX, &field[start ..]) {
                    let id = &field[start .. start + end];
                    if include_sets.contains(id) {
                        if exclude_aho.is_none() {
                            return true;
                        }
                    };
                } else {
                    return false;
                };
            } else if include_sets.contains(field) {
                if exclude_aho.is_none() {
                    return true;
                };
            } else {
                return false;
            }
        } else if field_index == 3 {
            // Field 4 (LCA): remainder after the last tab
            field_start += tab_pos + 1;
            let lca;
            if let Some(pos) = memchr(b'\t', &line[field_start ..]) {
                lca = &line[field_start .. (field_start + pos)]
            } else {
                lca = &line[field_start ..]
            };
            if let Some(ref exclude_matcher) = exclude_aho {
                return exclude_matcher.find(lca).is_none();
            }
        }
        field_index += 1;
        field_start += tab_pos + 1;
    }
    false
}
