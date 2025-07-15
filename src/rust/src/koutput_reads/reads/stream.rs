use anyhow::{anyhow, Result};
use bytes::{BufMut, Bytes};
use crossbeam_channel::Sender;
use libdeflater::Compressor;
use memchr::memchr;
use rustc_hash::FxHashMap as HashMap;

use crate::utils::*;

pub(in crate::koutput_reads::reads) struct KoutreadStream<H> {
    sender: Sender<Vec<u8>>,
    buffer: Vec<u8>,
    chunk_bytes: usize,
    tags: HashMap<Bytes, Bytes>,
    compressor: Option<Compressor>,
    handler: H,
}

pub(in crate::koutput_reads::reads) trait RecordHandler {
    type Record;

    /// Validate record correctness using external length info
    fn setup(&mut self, length: &Bytes, record: &Self::Record) -> Result<()>;
    fn seq_len(&self, record: &Self::Record) -> usize;
    fn qual_len(&self, record: &Self::Record) -> usize;
    /// Extract tags from the record (description and/or sequence)
    fn write_tags(&self, tags: &mut HashMap<Bytes, Bytes>, record: &Self::Record) -> Result<()>;

    /// Write sequence part to buffer
    fn write_seq(&self, buf: &mut Vec<u8>, record: &Self::Record);

    /// Write quality part to buffer
    fn write_qual(&self, buf: &mut Vec<u8>, record: &Self::Record);
}

impl<H> KoutreadStream<H>
where
    H: RecordHandler,
{
    #[allow(dead_code)]
    pub(in crate::koutput_reads::reads) fn new(sender: Sender<Vec<u8>>, handler: H) -> Self {
        Self::with_capacity(BLOCK_SIZE, sender, handler)
    }

    pub(in crate::koutput_reads::reads) fn with_capacity(
        capacity: usize,
        sender: Sender<Vec<u8>>,
        handler: H,
    ) -> Self {
        Self {
            sender,
            buffer: Vec::with_capacity(capacity),
            chunk_bytes: capacity,
            tags: HashMap::with_capacity_and_hasher(2, rustc_hash::FxBuildHasher),
            compressor: None,
            handler,
        }
    }

    pub(in crate::koutput_reads::reads) fn set_compressor(
        &mut self,
        compressor: Option<Compressor>,
    ) {
        self.compressor = compressor;
    }

    pub(in crate::koutput_reads::reads) fn process_record(
        &mut self,
        taxid: &Bytes,
        lca: &Bytes,
        length: &Bytes,
        record: &H::Record,
    ) -> Result<()> {
        self.handler.setup(length, record)?;

        // Extract tags from description field if any
        self.handler.write_tags(&mut self.tags, record)?;

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
            + self.handler.seq_len(record)
            + self.handler.qual_len(record)
            + 5;

        // If not enough buffer space, flush
        if self.buffer.capacity() - self.buffer.len() < len {
            let mut pack = Vec::with_capacity(self.chunk_bytes);
            std::mem::swap(&mut self.buffer, &mut pack);
            self.send(pack)?;
        }

        #[cfg(debug_assertions)]
        let start = self.buffer.len();

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
        self.handler.write_seq(&mut self.buffer, record);
        self.buffer.put_u8(b'\t');
        self.handler.write_qual(&mut self.buffer, record);
        self.buffer.put_u8(b'\n');

        // Debug assertion to verify length matches expectation
        #[cfg(debug_assertions)]
        assert_eq!(
            self.buffer.len() - start,
            len,
            "Buffer write length mismatch: expected {}, actual {}",
            len,
            self.buffer.len() - start
        );

        self.tags.clear();

        Ok(())
    }

    pub(in crate::koutput_reads::reads) fn send(&mut self, mut pack: Vec<u8>) -> Result<()> {
        // Compress if gzip file
        if let Some(compressor) = &mut self.compressor {
            pack = gzip_pack(&pack, compressor)?
        }

        // Send compressed or raw bytes to writer
        self.sender.send(pack).map_err(|e| {
            anyhow!(
                "(Parser) Failed to send parsed lines to Writer thread: {}",
                e
            )
        })?;
        Ok(())
    }

    pub(in crate::koutput_reads::reads) fn flush_buffer(mut self) -> Result<()> {
        let pack = std::mem::take(&mut self.buffer);
        self.send(pack)
    }
}

pub(in crate::koutput_reads::reads) fn extract_tags_from_desc(
    tags: &mut HashMap<Bytes, Bytes>,
    desc: &Option<Bytes>,
) {
    if let Some(desc) = desc {
        if let Some(start) = TAG_PREFIX_FINDER.find(desc) {
            if let Some(end) = memchr(TAG_SUFFIX, &desc[start ..]) {
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
                            tags.insert(tag, sequence);
                        }
                    }
                }
            }
        }
    }
}
