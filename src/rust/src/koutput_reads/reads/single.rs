use std::io::Write;
use std::path::Path;

use anyhow::{anyhow, Result};
use bytes::BufMut;
use bytes::Bytes;
use bytes::BytesMut;
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::Compressor;
use memchr::memchr;
use rustc_hash::FxHashMap as HashMap;

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
    tag_ranges: Option<&TagRanges>,
    chunk_size: usize,
    buffer_size: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let input: &Path = input_path.as_ref();
    let output: &Path = output_path.as_ref();
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
            let mut writer = new_writer(output, buffer_size, output_bar)?;

            // Iterate over each received batch of records
            for chunk in writer_rx {
                writer.write_all(&chunk).map_err(|e| {
                    anyhow!("(Writer) Failed to write FastqRecord to output: {}", e)
                })?;
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let mut parser_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                let mut stream = KoutreadStream::new(tx);
                stream.set_tag_ranges(tag_ranges);
                while let Ok(records) = rx.recv() {
                    for record in records {
                        // sequence length, taxid, lca
                        if let Some((length, taxid, lca)) = koutmap.get(&record.id) {
                            if memchr(b':', length).is_some() {
                                return Err(anyhow!(
                                    "Invalid input: paired-end format detected in kraken2 output, but only single-end reads were provided"
                                ));
                            }
                            let expected_len = std::str::from_utf8(length)?
                                .parse::<usize>()
                                .map_err(|e| {
                                    anyhow!(
                                        "Invalid length in koutput for ID {:?}: {}",
                                        record.id, e
                                    )
                                })?;
                            if expected_len != record.seq.len() {
                                return Err(anyhow!(
                                    "Invalid input: expected sequence length in kraken2 output {}, got {}",
                                    expected_len,
                                    record.seq.len()
                                ));
                            }
                            stream.process_record(taxid, lca, &record)?;
                        }
                    }
                }
                stream.flush_buffer().map_err(|e| {
                    anyhow!("(Parser) Failed to flush parsed lines to Writer thread: {}", e)
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
                FastqReader::with_capacity(buffer_size, new_reader(input, buffer_size, input_bar)?);
            let mut reader_tx = BatchSender::with_capacity(chunk_size, reader_tx);
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

struct KoutreadStream<'a> {
    sender: Sender<Vec<u8>>,
    buffer: Vec<u8>,
    tags: HashMap<Bytes, Bytes>,
    tag_ranges: Option<&'a TagRanges>,
    compressor: Option<Compressor>,
}

impl<'a> KoutreadStream<'a> {
    fn new(sender: Sender<Vec<u8>>) -> Self {
        Self {
            sender,
            buffer: Vec::with_capacity(BLOCK_SIZE),
            tags: HashMap::with_capacity_and_hasher(2, rustc_hash::FxBuildHasher),
            tag_ranges: None,
            compressor: None,
        }
    }

    fn set_tag_ranges(&mut self, tag_ranges: Option<&'a TagRanges>) {
        self.tag_ranges = tag_ranges;
    }

    fn set_compressor(&mut self, compressor: Option<Compressor>) {
        self.compressor = compressor;
    }

    fn process_record(
        &mut self,
        taxid: &Bytes,
        lca: &Bytes,
        record: &FastqRecord<Bytes>,
    ) -> Result<()> {
        // Extract tags from description field if any
        self.extract_tags_from_desc(&record.desc);

        // Extract tags from current sequence
        self.extract_tags_from_seq(&record.seq)?;

        // Precompute required space: taxid + tags + lca + seq + qual + 4 tabs + 1 newline
        let len = taxid.len()
                + self
                    .tags
                    .iter()
                    // tag:sequence
                    .map(|(tag, sequence)| tag.len() + 1 + sequence.len())
                    .sum::<usize>()
                + self.tags.len().saturating_sub(1) // space separators if more than 1
            + lca.len()
            + record.seq.len()
            + record.qual.len()
            + 5;

        // If not enough buffer space, flush
        if self.buffer.capacity() - self.buffer.len() < len {
            let mut pack = Vec::with_capacity(BLOCK_SIZE);
            std::mem::swap(&mut self.buffer, &mut pack);
            self.send(pack)?;
        }

        // Write fields to buffer: taxid \t tags \t lca \t seq \t qual \n
        self.buffer.extend_from_slice(taxid);
        self.buffer.put_u8(b'\t');

        for (i, (tag, sequence)) in self.tags.iter().enumerate() {
            if i > 0 {
                self.buffer.put_u8(b' ');
            }
            self.buffer.extend_from_slice(&tag);
            self.buffer.put_u8(b':');
            self.buffer.extend_from_slice(sequence);
        }
        self.buffer.put_u8(b'\t');
        self.buffer.extend_from_slice(lca);
        self.buffer.put_u8(b'\t');
        self.buffer.extend_from_slice(&record.seq);
        self.buffer.put_u8(b'\t');
        self.buffer.extend_from_slice(&record.qual);
        self.buffer.put_u8(b'\n');
        self.tags.clear();

        Ok(())
    }

    fn send(&mut self, mut pack: Vec<u8>) -> Result<()> {
        // Compress if gzip file
        if let Some(compressor) = &mut self.compressor {
            pack = gzip_pack(&pack, compressor)?
        }

        // Send compressed or raw bytes to writer
        self.sender.send(pack)?;
        Ok(())
    }

    fn flush_buffer(mut self) -> Result<()> {
        let pack = std::mem::take(&mut self.buffer);
        self.send(pack)
    }

    fn extract_tags_from_desc(&mut self, desc: &Option<Bytes>) {
        if let Some(desc) = desc {
            if let Some(start) = TAG_PREFIX_FINDER.find(desc) {
                if let Some(end) = memchr(b'}', &desc[start ..]) {
                    // Inside `RSAHMI{}`
                    let buf = &desc[start + TAG_PREFIX.len() .. start + end];

                    // Parse as tag:seq:tag:seq:...
                    let mut pos = 0;
                    while pos < buf.len() {
                        if let Some(separator) = memchr(b':', &buf[pos ..]) {
                            let start = pos; // field start
                            pos += separator + 1;
                            if pos < buf.len() {
                                let tag = desc.slice_ref(&buf[start .. start + separator]);
                                let sequence;
                                if let Some(end) = memchr(b':', &buf[pos ..]) {
                                    sequence = desc.slice_ref(&buf[pos .. pos + end]);
                                    pos += end + 1;
                                } else {
                                    sequence = desc.slice_ref(&buf[pos .. buf.len()]);
                                    pos = buf.len();
                                }
                                self.tags.insert(tag, sequence);
                            }
                        }
                    }
                }
            }
        }
    }

    fn extract_tags_from_seq(&mut self, seq: &Bytes) -> Result<()> {
        if let Some(tag_ranges) = self.tag_ranges {
            let tag_map = tag_ranges.map_sequences(seq).map_err(|e| {
                anyhow!(
                    "(Parser) Failed to send parsed lines to Writer thread: {}",
                    e
                )
            })?;
            for (tag, sequences) in tag_map {
                let mut seq = BytesMut::with_capacity(sequences.iter().map(|s| s.len()).sum());
                for sequence in sequences {
                    seq.extend_from_slice(sequence);
                }
                self.tags.insert(tag, seq.freeze());
            }
        }
        Ok(())
    }
}
