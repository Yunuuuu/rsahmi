use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, Result};
use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::{CompressionLvl, Compressor};
use rustc_hash::FxHashMap as HashMap;

use super::stream::extract_tags_from_desc;
use super::stream::RecordHandler;
use crate::batchsender::BatchSender;
use crate::fastq_reader::*;
use crate::parser::fastq::FastqRecord;
use crate::seq_tag::*;
use crate::utils::*;

pub(crate) fn parse_single_read<P: AsRef<Path> + ?Sized>(
    koutmap: &HashMap<Bytes, (Bytes, Bytes, Bytes)>,
    input_path: &P,
    input_bar: Option<ProgressBar>,
    output_path: &P,
    output_bar: Option<ProgressBar>,
    tag_ranges: &Option<TagRanges>,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: CompressionLvl,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let input: &Path = input_path.as_ref();
    let output: &Path = output_path.as_ref();
    let gzip = gz_compressed(output);
    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (writer_tx, writer_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = new_channel(nqueue);

        let (reader_tx, reader_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            if let Some(bar) = &output_bar {
                bar.tick();
            }
            let mut writer = BufWriter::with_capacity(chunk_bytes, new_writer(output, output_bar)?);

            // Iterate over each received batch of records
            for chunk in writer_rx {
                writer
                    .write_all(&chunk)
                    .map_err(|e| anyhow!("(Writer) Failed to write to output: {}", e))?;
            }
            writer
                .flush()
                .map_err(|e| anyhow!("(Writer) Failed to flush writer: {}", e))?;
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let mut parser_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                let record_handler = SinlgeRecordHandle::new(tag_ranges);
                let mut stream = crate::koutput_reads::reads::stream::KoutreadStream::with_capacity(
                    chunk_bytes,
                    tx,
                    record_handler,
                );
                if gzip {
                    let compressor = Compressor::new(compression_level);
                    stream.set_compressor(Some(compressor));
                }
                while let Ok(records) = rx.recv() {
                    for record in records {
                        if let Some((length, taxid, lca)) = koutmap.get(&record.id) {
                            stream.process_record(taxid, lca, length, &record)?;
                        }
                    }
                }
                stream.flush_buffer().map_err(|e| {
                    anyhow!(
                        "(Parser) Failed to flush parsed lines to Writer thread: {}",
                        e
                    )
                })?;
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
                .map_err(|e| anyhow!("(Reader) Error while reading FASTQ record: {}", e))?
            {
                reader_tx.send(record).map_err(|e| {
                    anyhow!(
                        "(Reader) Failed to send FASTQ record to Parser thread: {}",
                        e
                    )
                })?;
            }
            reader_tx
                .flush()
                .map_err(|e| anyhow!("(Reader) Failed to flush records to Parser thread: {}", e))?;
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

struct SinlgeRecordHandle<'a> {
    tag_ranges: &'a Option<TagRanges>,
}

impl<'a> SinlgeRecordHandle<'a> {
    fn new(tag_ranges: &'a Option<TagRanges>) -> Self {
        Self { tag_ranges }
    }
}

impl<'a> RecordHandler for SinlgeRecordHandle<'a> {
    type Record = FastqRecord<Bytes>;
    fn setup(&mut self, length: &Bytes, record: &Self::Record) -> Result<()> {
        if memchr::memchr(b':', length).is_some() {
            return Err(anyhow!(
                "Invalid input: paired-end format detected in kraken2 output, but only single-end reads were provided"
            ));
        }

        let expected_len = std::str::from_utf8(length)?
            .trim()
            .parse::<usize>()
            .map_err(|e| anyhow!("Invalid length in koutput for ID {:?}: {}", record.id, e))?;

        if expected_len != record.seq.len() {
            return Err(anyhow!(
                "Sequence length mismatch: expected {}, got {}",
                expected_len,
                record.seq.len()
            ));
        }

        Ok(())
    }
    fn seq_len(&self, record: &Self::Record) -> usize {
        record.seq.len()
    }

    fn qual_len(&self, record: &Self::Record) -> usize {
        record.qual.len()
    }

    fn write_tags(&self, tags: &mut HashMap<Bytes, Bytes>, record: &Self::Record) -> Result<()> {
        // 1. Extract from description
        extract_tags_from_desc(tags, &record.desc);

        // 2. Extract from sequence using tag_ranges
        if let Some(tag_ranges) = self.tag_ranges {
            let tag_map = tag_ranges.map_sequences(&record.seq).map_err(|e| {
                anyhow!("(TagExtractor) Failed to extract tags from sequence: {}", e)
            })?;
            for (tag, sequences) in tag_map {
                let mut seq = BytesMut::with_capacity(sequences.iter().map(|s| s.len()).sum());
                for sequence in sequences {
                    seq.extend_from_slice(sequence);
                }
                tags.insert(tag, seq.freeze());
            }
        }

        Ok(())
    }

    fn write_seq(&self, buf: &mut Vec<u8>, record: &Self::Record) {
        buf.extend_from_slice(&record.seq);
    }

    fn write_qual(&self, buf: &mut Vec<u8>, record: &Self::Record) {
        buf.extend_from_slice(&record.qual);
    }
}
