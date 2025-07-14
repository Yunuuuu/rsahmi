use std::io::Write;
use std::iter::zip;
use std::path::Path;

use anyhow::{anyhow, Result};
use bytes::BufMut;
use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::{CompressionLvl, Compressor};
use rustc_hash::FxHashMap as HashMap;

use super::stream::extract_tags_from_desc;
use super::stream::RecordHandler;
use crate::batchsender::BatchSender;
use crate::fastq_reader::*;
use crate::koutput_reads::reads::stream::KoutreadStream;
use crate::parser::fastq::FastqRecord;
use crate::seq_tag::*;
use crate::utils::*;

pub(crate) fn parse_paired_read<P: AsRef<Path> + ?Sized>(
    koutmap: &HashMap<Bytes, (Bytes, Bytes, Bytes)>,
    input1_path: &P,
    input1_bar: Option<ProgressBar>,
    input2_path: &P,
    input2_bar: Option<ProgressBar>,
    output_path: &P,
    output_bar: Option<ProgressBar>,
    tag_ranges1: &Option<TagRanges>,
    tag_ranges2: &Option<TagRanges>,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: CompressionLvl,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let output: &Path = output_path.as_ref();
    let gzip = gz_compressed(output);
    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (writer_tx, writer_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = new_channel(nqueue);

        let (reader_tx, reader_rx): (
            Sender<(Vec<FastqRecord<Bytes>>, Vec<FastqRecord<Bytes>>)>,
            Receiver<(Vec<FastqRecord<Bytes>>, Vec<FastqRecord<Bytes>>)>,
        ) = new_channel(nqueue);
        let (reader1_tx, reader1_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = new_channel(nqueue);
        let (reader2_tx, reader2_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            let mut writer = new_writer(output, output_bar)?;

            // Iterate over each received batch of records
            for chunk in writer_rx {
                writer
                    .write_all(&chunk)
                    .map_err(|e| anyhow!("(Writer) Failed to write to output: {}", e))?;
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let mut parser_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                let record_handler = PairedRecordHandle::new(tag_ranges1, tag_ranges2);
                let mut stream = KoutreadStream::with_capacity(chunk_bytes, tx, record_handler);
                if gzip {
                    let compressor = Compressor::new(compression_level);
                    stream.set_compressor(Some(compressor));
                }
                while let Ok((records1, records2)) = rx.recv() {
                    // Initialize a thread-local batch sender for matching records
                    for (record1, record2) in zip(records1, records2) {
                        if record1.id != record2.id {
                            return Err(anyhow!(
                                "FASTQ pairing error: record1 ID = {}, record2 ID = {}. These records do not match and cannot be paired.",
                                String::from_utf8_lossy(&record1.id),
                                String::from_utf8_lossy(&record2.id)
                            ));
                        }
                        if let Some((length, taxid, lca)) = koutmap.get(&record1.id) {
                            stream.process_record(taxid, lca, length, &(record1, record2))?;
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
            loop {
                let (records1, records2) = match (reader1_rx.recv(), reader2_rx.recv()) {
                    (Ok(rec1), Ok(rec2)) => (rec1, rec2),
                    (Err(_), Ok(_)) => {
                        return Err(anyhow!(
                            "(Reader collect) FASTQ pairing error: read1 channel closed before read2"
                        ));
                    }
                    (Ok(_), Err(_)) => {
                        return Err(anyhow!(
                            "(Reader collect) FASTQ pairing error: read2 channel closed before read1"
                        ));
                    }
                    (Err(_), Err(_)) => {
                        break;
                    }
                };
                if records1.len() != records2.len() {
                    return Err(anyhow!("(Reader collect) FASTQ pairing error: record count mismatch (read1: {}, read2: {})", records1.len(), records2.len()));
                }
                reader_tx.send((records1, records2)).map_err(|e| {
                    anyhow!(
                        "(Reader collect) Failed to send send parsed record pair to Parser thread: {}",
                        e
                    )
                })?;
            }
            Ok(())
        });

        let input1: &Path = input1_path.as_ref();
        let reader1_handle = scope.spawn(move || -> Result<()> {
            let mut reader = FastqReader::with_capacity(
                BUFFER_SIZE,
                new_reader(input1, BUFFER_SIZE, input1_bar)?,
            );
            let mut thread_tx = BatchSender::with_capacity(batch_size, reader1_tx);
            while let Some(record) = reader
                .read_record()
                .map_err(|e| anyhow!("(Reader1) Error while reading FASTQ record: {}", e))?
            {
                thread_tx.send(record).map_err(|e| {
                    anyhow!(
                        "(Reader1) Failed to send FASTQ record to reader collect thread: {}",
                        e
                    )
                })?;
            }
            thread_tx.flush().map_err(|e| {
                anyhow!(
                    "(Reader1) Failed to flush records to reader collect thread: {}",
                    e
                )
            })?;
            Ok(())
        });

        let input2: &Path = input2_path.as_ref();
        let reader2_handle = scope.spawn(move || -> Result<()> {
            let mut reader = FastqReader::with_capacity(
                BUFFER_SIZE,
                new_reader(input2, BUFFER_SIZE, input2_bar)?,
            );
            let mut thread_tx = BatchSender::with_capacity(batch_size, reader2_tx);
            while let Some(record) = reader
                .read_record()
                .map_err(|e| anyhow!("(Reader2) Error while reading FASTQ record: {}", e))?
            {
                thread_tx.send(record).map_err(|e| {
                    anyhow!(
                        "(Reader2) Failed to send FASTQ record to reader collect thread: {}",
                        e
                    )
                })?;
            }
            thread_tx.flush().map_err(|e| {
                anyhow!(
                    "(Reader2) Failed to flush records to reader collect thread: {}",
                    e
                )
            })?;
            Ok(())
        });

        // ─── Join Threads and Propagate Errors ────────────────
        writer_handle
            .join()
            .map_err(|e| anyhow!("(Writer dispatch) thread panicked: {:?}", e))??;

        for handler in parser_handles {
            handler
                .join()
                .map_err(|e| anyhow!("(Parser) thread panicked: {:?}", e))??;
        }
        reader_handle
            .join()
            .map_err(|e| anyhow!("(Reader collect) thread panicked: {:?}", e))??;
        reader1_handle
            .join()
            .map_err(|e| anyhow!("(Reader1) thread panicked: {:?}", e))??;
        reader2_handle
            .join()
            .map_err(|e| anyhow!("(Reader2) thread panicked: {:?}", e))??;
        Ok(())
    })
}

struct PairedRecordHandle<'a> {
    tag_ranges1: &'a Option<TagRanges>,
    tag_ranges2: &'a Option<TagRanges>,
    pair: bool,
}

impl<'a> PairedRecordHandle<'a> {
    fn new(tag_ranges1: &'a Option<TagRanges>, tag_ranges2: &'a Option<TagRanges>) -> Self {
        Self {
            tag_ranges1,
            tag_ranges2,
            pair: false,
        }
    }
}

impl<'a> RecordHandler for PairedRecordHandle<'a> {
    type Record = (FastqRecord<Bytes>, FastqRecord<Bytes>);
    fn setup(&mut self, length: &Bytes, record: &Self::Record) -> Result<()> {
        if let Some(separator) = memchr::memchr(b':', length) {
            let (l1, l2) = (&length[.. separator], &length[separator + 1 ..]);
            let len1 = std::str::from_utf8(l1)?
                .trim()
                .parse::<usize>()
                .map_err(|e| {
                    anyhow!(
                        "(Paired read1) Invalid length in koutput for ID {:?}: {}",
                        record.0.id,
                        e
                    )
                })?;
            let len2 = std::str::from_utf8(l2)?
                .trim()
                .parse::<usize>()
                .map_err(|e| {
                    anyhow!(
                        "(Paired Read2) Invalid length in koutput for ID {:?}: {}",
                        record.1.id,
                        e
                    )
                })?;
            if record.0.seq.len() != len1 || record.1.seq.len() != len2 {
                return Err(anyhow!(
                    "(Read1 and Read2) Invalid input: expected lengths {}:{}, got {}:{}",
                    len1,
                    len2,
                    record.0.seq.len(),
                    record.1.seq.len()
                ));
            }
            self.pair = true;
        } else {
            let expected_len = std::str::from_utf8(length)?
                .trim()
                .parse::<usize>()
                .map_err(|e| {
                    anyhow!(
                        "(Read2) Invalid length in koutput for ID {:?}: {}",
                        record.1.id,
                        e
                    )
                })?;
            if expected_len != record.1.seq.len() {
                return Err(anyhow!(
                    "(Read2) Sequence length mismatch: expected {}, got {}",
                    expected_len,
                    record.1.seq.len()
                ));
            }
            self.pair = false;
        };

        Ok(())
    }

    fn seq_len(&self, record: &Self::Record) -> usize {
        if self.pair {
            record.0.seq.len() + 1 + record.1.seq.len()
        } else {
            record.1.seq.len()
        }
    }

    fn qual_len(&self, record: &Self::Record) -> usize {
        if self.pair {
            record.0.qual.len() + 1 + record.1.qual.len()
        } else {
            record.1.qual.len()
        }
    }

    fn write_tags(&self, tags: &mut HashMap<Bytes, Bytes>, record: &Self::Record) -> Result<()> {
        // 1. Extract from description
        extract_tags_from_desc(tags, &record.0.desc);
        extract_tags_from_desc(tags, &record.1.desc);

        // 2. Extract from sequence using tag_ranges
        let tag_map = match (&self.tag_ranges1, &self.tag_ranges2) {
            (Some(tag_ranges), None) => tag_ranges.map_sequences(&record.0.seq).map_err(|e| {
                anyhow!("(TagExtractor) Failed to extract tags from sequence: {}", e)
            })?,
            (None, Some(tag_ranges)) => tag_ranges.map_sequences(&record.1.seq).map_err(|e| {
                anyhow!("(TagExtractor) Failed to extract tags from sequence: {}", e)
            })?,
            (Some(tag_ranges1), Some(tag_ranges2)) => {
                let mut tag_map = tag_ranges1.map_sequences(&record.0.seq).map_err(|e| {
                    anyhow!("(TagExtractor) Failed to extract tags from sequence: {}", e)
                })?;
                let tag_map2 = tag_ranges2.map_sequences(&record.1.seq).map_err(|e| {
                    anyhow!("(TagExtractor) Failed to extract tags from sequence: {}", e)
                })?;

                // Merge tag→sequence entries
                for (tag, sequence) in tag_map2 {
                    if let Some(v) = tag_map.get_mut(&tag) {
                        v.extend(sequence); // read1 first, read2 second
                    } else {
                        tag_map.insert(tag, sequence);
                    }
                }
                tag_map
            }
            (None, None) => return Ok(()),
        };
        for (tag, sequences) in tag_map {
            let mut seq = BytesMut::with_capacity(sequences.iter().map(|s| s.len()).sum());
            for sequence in sequences {
                seq.extend_from_slice(sequence);
            }
            tags.insert(tag, seq.freeze());
        }
        Ok(())
    }

    fn write_seq(&self, buf: &mut Vec<u8>, record: &Self::Record) {
        if self.pair {
            buf.extend_from_slice(&record.0.seq);
            buf.put_u8(b'\t');
        }
        buf.extend_from_slice(&record.1.seq);
    }

    fn write_qual(&self, buf: &mut Vec<u8>, record: &Self::Record) {
        if self.pair {
            buf.extend_from_slice(&record.0.qual);
            buf.put_u8(b'\t');
        }
        buf.extend_from_slice(&record.1.qual);
    }
}
