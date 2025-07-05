use std::cell::Cell;
use std::io::Read;

use anyhow::{anyhow, Result};
use bytes::Bytes;

use crate::parser::fastq::{FastqContainer, FastqParseError, FastqRecord};
use crate::reader::bytes::{BytesBreaksReader, BytesReader};

pub(crate) struct FastqBytesChunkReader<R: Read> {
    line_offset: usize,
    chunk_size: usize,
    reader: BytesReader<R>,
}

impl<R: Read> FastqBytesChunkReader<R> {
    #[allow(dead_code)]
    pub(crate) fn new(reader: BytesReader<R>) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub(crate) fn with_capacity(capacity: usize, reader: BytesReader<R>) -> Self {
        Self {
            line_offset: 0,
            chunk_size: capacity,
            reader,
        }
    }

    pub(crate) fn chunk_reader(
        &mut self,
    ) -> Result<Option<FastqBytesReader<std::vec::IntoIter<usize>>>> {
        // read files
        let buf = match self.reader.read_bytes(self.chunk_size)? {
            Some(buf) => buf,
            None => return Ok(None),
        };

        // Find all newline character positions in the chunk
        let mut newlines = buf.newlines();

        // If the chunk contains fewer than 4 lines, the FASTQ record is incomplete
        // so we expand the buffer and try again with a larger chunk
        if newlines.len() < 4 && !buf.eof() {
            self.reader.take_leftover(buf.into_bytes()); // remainder for next round
                                                         // we always double the buffer size, if there no enough buffer to hold a whole fastq record
            self.chunk_size *= 2;
            return self.chunk_reader();
        }

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        if newlines.len() >= 4 {
            newlines.truncate(newlines.len() - (newlines.len() % 4));
        }

        // SAFETY: above code will ensure newlines won't be empty
        let endpoint = unsafe { *newlines.last().unwrap_unchecked() };
        let mut bytes = buf.into_bytes();
        let chunk;
        if endpoint < bytes.len() {
            chunk = bytes.split_to(endpoint + 1); // include `\n`
            self.reader.take_leftover(bytes);
        } else {
            chunk = bytes;
        }

        let line_offset = self.line_offset;
        self.line_offset += newlines.len(); // Increment the number of chunk lines
        let mut source = FastqBytesSource::new(chunk.freeze(), newlines.into_iter());
        source.set_offset(line_offset);
        if let Some(label) = self.reader.label() {
            source.set_label(label);
        }
        return Ok(Some(FastqBytesReader::new(source)));
    }
}

pub(crate) struct FastqBytesChunkPairedReader<R1: Read, R2: Read> {
    line_offset1: usize,
    line_offset2: usize,
    chunk_size: usize,
    reader1: BytesReader<R1>,
    reader2: BytesReader<R2>,
}

impl<R1: Read, R2: Read> FastqBytesChunkPairedReader<R1, R2> {
    #[allow(dead_code)]
    pub(crate) fn new(reader1: BytesReader<R1>, reader2: BytesReader<R2>) -> Self {
        Self::with_capacity(8 * 1024, reader1, reader2)
    }

    pub(crate) fn with_capacity(
        capacity: usize,
        reader1: BytesReader<R1>,
        reader2: BytesReader<R2>,
    ) -> Self {
        Self {
            line_offset1: 0,
            line_offset2: 0,
            chunk_size: capacity,
            reader1,
            reader2,
        }
    }

    pub(crate) fn chunk_reader(
        &mut self,
    ) -> Result<Option<BytesFastqPairedChunkReader<std::vec::IntoIter<usize>>>> {
        // --- Read data for read1 and append to existing leftover ---
        let buf1 = self.reader1.read_bytes(self.chunk_size)?;
        let buf2 = self.reader2.read_bytes(self.chunk_size)?;

        // --- EOF + desync checks ---
        let (buf1, buf2) = match (buf1, buf2) {
            (None, None) => return Ok(None),

            // read1 ended, check if read2 has only trailing whitespace
            (None, Some(ref bytes)) => {
                if bytes.as_slice().iter().all(|b| b.is_ascii_whitespace()) {
                    // read and discard trailing whitespace from stream
                    for byte_result in self.reader2.borrow_reader().bytes() {
                        let byte = byte_result?;
                        if !byte.is_ascii_whitespace() {
                            break;
                        }
                    }
                    return Ok(None);
                }
                // Possible true sync failure
                return Err(anyhow!("{}", FastqParseError::OutOfSync {
                    eof_label: self.reader1.label(),
                    continueal_label: self.reader2.label(),
                    eof_pos: self.line_offset1,
                    continueal_pos: self.line_offset2,
                }));
            }
            (Some(ref bytes), None) => {
                if bytes.as_slice().iter().all(|b| b.is_ascii_whitespace()) {
                    for byte_result in self.reader1.borrow_reader().bytes() {
                        let byte = byte_result?;
                        if !byte.is_ascii_whitespace() {
                            break;
                        }
                    }
                    return Ok(None);
                }
                return Err(anyhow!("{}", FastqParseError::OutOfSync {
                    eof_label: self.reader2.label(),
                    continueal_label: self.reader1.label(),
                    eof_pos: self.line_offset2,
                    continueal_pos: self.line_offset1,
                }));
            }
            (Some(bytes1), Some(bytes2)) => (bytes1, bytes2),
        };

        // --- Find newlines in both chunks ---
        // read1
        let mut newlines1 = buf1.newlines();

        // read2
        let mut newlines2 = buf2.newlines();

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        // --- Determine how many records we can safely return (must be multiple of 4 lines) ---
        let count = usize::min(newlines1.len(), newlines2.len());
        if count < 4 {
            if !buf1.eof() && !buf2.eof() {
                // Not enough lines yet, double chunk size for next round
                self.reader1.take_leftover(buf1.into_bytes());
                self.reader2.take_leftover(buf2.into_bytes());
                self.chunk_size *= 2;
                return self.chunk_reader();
            } else {
                // Partial record at EOF — invalid
                // If one file only contain few lines, it means it has incompleted record
                let (label, pos, record) = if buf1.eof() {
                    (
                        self.reader1.label(),
                        self.line_offset1,
                        String::from_utf8_lossy(buf1.as_slice()).to_string(),
                    )
                } else {
                    (
                        self.reader2.label(),
                        self.line_offset2,
                        String::from_utf8_lossy(buf2.as_slice()).to_string(),
                    )
                };
                return Err(anyhow!("{}", FastqParseError::IncompleteRecord {
                    label,
                    record,
                    pos
                }));
            }
        }

        // --- Cut to exact multiple of 4 ---
        let nkeep = count - (count % 4);
        newlines1.truncate(nkeep);
        newlines2.truncate(nkeep);

        // --- Split buf1 using endpoint newline ---
        let endpoint1 = unsafe { *newlines1.last().unwrap_unchecked() };
        let mut bytes = buf1.into_bytes();
        let chunk1;
        if endpoint1 < bytes.len() {
            chunk1 = bytes.split_to(endpoint1 + 1); // include `\n`
            self.reader1.take_leftover(bytes);
        } else {
            chunk1 = bytes;
        }

        // --- Split buf2 likewise ---
        let endpoint2 = unsafe { *newlines2.last().unwrap_unchecked() };
        let mut bytes = buf2.into_bytes();
        let chunk2;
        if endpoint2 < bytes.len() {
            chunk2 = bytes.split_to(endpoint2 + 1); // include `\n`
            self.reader2.take_leftover(bytes);
        } else {
            chunk2 = bytes;
        }

        // --- Update line offset bookkeeping ---
        let line_offset1 = self.line_offset1;
        let line_offset2 = self.line_offset2;
        self.line_offset1 += newlines1.len();
        self.line_offset2 += newlines2.len();

        // --- Wrap into sources with offsets and labels ---
        let mut source1 = FastqBytesSource::new(chunk1.freeze(), newlines1.into_iter());
        source1.set_offset(line_offset1);
        if let Some(label1) = self.reader1.label() {
            source1.set_label(label1);
        }

        let mut source2 = FastqBytesSource::new(chunk2.freeze(), newlines2.into_iter());
        source2.set_offset(line_offset2);
        if let Some(label2) = self.reader2.label() {
            source2.set_label(label2);
        }
        Ok(Some(BytesFastqPairedChunkReader::new(
            FastqBytesReader::new(source1),
            FastqBytesReader::new(source2),
        )))
    }
}

#[derive(Debug)]
pub(crate) struct FastqBytesSource<I: Iterator<Item = usize>> {
    // we need hold long-lived references to the chunk, so we use
    // normal reference, and wrap `offset` into Cell
    // so we can change them
    offset: Cell<usize>,
    reader: BytesBreaksReader<I>,
}

impl<I: Iterator<Item = usize>> FastqBytesSource<I> {
    fn new(bytes: Bytes, breaks: I) -> Self {
        Self {
            offset: Cell::new(0),
            reader: BytesBreaksReader::with_breaks(breaks, bytes),
        }
    }

    fn read_line(&self) -> Option<&[u8]> {
        let out = self.reader.read_line();
        if out.is_some() {
            self.set_offset(self.offset.get() + 1);
        }
        out
    }

    fn label(&self) -> Option<&'static str> {
        self.reader.label()
    }

    fn set_label(&mut self, label: &'static str) {
        self.reader.set_label(label)
    }

    #[allow(dead_code)]
    fn unset_label(&mut self) {
        self.reader.unset_label()
    }

    // when split files into multiple chunks, we need adjust position manually.
    fn set_offset(&self, offset: usize) {
        self.offset.set(offset)
    }

    fn get_offset(&self) -> usize {
        self.offset.get()
    }

    fn borrow_bytes(&self) -> &Bytes {
        self.reader.borrow_bytes()
    }
}

#[derive(Debug)]
pub(crate) struct FastqBytesReader<I: Iterator<Item = usize>> {
    source: FastqBytesSource<I>,
}

impl<I: Iterator<Item = usize>> FastqBytesReader<I> {
    pub(crate) fn new(source: FastqBytesSource<I>) -> Self {
        Self { source }
    }

    pub(crate) fn read_record<'inner, 'outer>(
        &'outer self,
        container: &mut FastqContainer<'inner>,
    ) -> Result<Option<FastqRecord<Bytes>>, FastqParseError>
    where
        'outer: 'inner,
    {
        // Try reading the 1st line (head). If EOF, return None.
        if self.read_head(container)?.is_none() {
            return Ok(None);
        };
        self.read_tail(container)?;
        Ok(Some(self.build(container)))
    }

    fn read_head<'inner, 'outer>(
        &'outer self,
        container: &mut FastqContainer<'inner>,
    ) -> Result<Option<()>, FastqParseError>
    where
        'outer: 'inner,
    {
        // Try reading the 1st line (head). If EOF, return None.
        loop {
            match self.source.read_line() {
                Some(line) => {
                    // Skip empty or all-whitespace lines
                    if line.is_empty() || line.iter().all(|b| b.is_ascii_whitespace()) {
                        continue;
                    }
                    container.parse_head(
                        line,
                        self.source.label(),
                        self.source.get_offset() - 1,
                    )?;
                    return Ok(Some(()));
                }
                None => return Ok(None), // EOF
            }
        }
    }

    fn read_tail<'inner, 'outer>(
        &'outer self,
        container: &mut FastqContainer<'inner>,
    ) -> Result<(), FastqParseError>
    where
        'outer: 'inner,
    {
        // 2nd line (sequence) must exist. Otherwise, incomplete record.
        if let Some(line) = self.source.read_line() {
            container.parse_seq(line, self.source.label(), self.source.get_offset() - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label(),
                record: container.to_string(),
                pos: self.source.get_offset(),
            });
        }

        // 3rd line (separator). Must exist.
        if let Some(line) = self.source.read_line() {
            container.parse_sep(line, self.source.label(), self.source.get_offset() - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label(),
                record: container.to_string(),
                pos: self.source.get_offset(),
            });
        }

        // 4th line (quality). Must exist.
        if let Some(line) = self.source.read_line() {
            container.parse_qual(line, self.source.label(), self.source.get_offset() - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label(),
                record: container.to_string(),
                pos: self.source.get_offset(),
            });
        }
        Ok(())
    }

    fn build<'inner, 'outer>(
        &'outer self,
        container: &mut FastqContainer<'inner>,
    ) -> FastqRecord<Bytes>
    where
        'outer: 'inner,
    {
        let out = unsafe {
            FastqRecord::new(
                self.source
                    .borrow_bytes()
                    .slice_ref(container.id().unwrap_unchecked()),
                container
                    .desc()
                    .unwrap_unchecked()
                    .map(|slice| self.source.borrow_bytes().slice_ref(slice)),
                self.source
                    .borrow_bytes()
                    .slice_ref(container.seq().unwrap_unchecked()),
                self.source
                    .borrow_bytes()
                    .slice_ref(container.sep().unwrap_unchecked()),
                self.source
                    .borrow_bytes()
                    .slice_ref(container.qual().unwrap_unchecked()),
            )
        };
        container.reset();
        out
    }
}

pub(crate) struct BytesFastqPairedChunkReader<I: Iterator<Item = usize>> {
    reader1: FastqBytesReader<I>,
    reader2: FastqBytesReader<I>,
}

impl<I: Iterator<Item = usize>> BytesFastqPairedChunkReader<I> {
    pub(crate) fn new(reader1: FastqBytesReader<I>, reader2: FastqBytesReader<I>) -> Self {
        Self { reader1, reader2 }
    }

    pub(crate) fn read_record<'outer, 'inner1, 'inner2>(
        &'outer self,
        container1: &mut FastqContainer<'inner1>,
        container2: &mut FastqContainer<'inner2>,
    ) -> Result<Option<(FastqRecord<Bytes>, FastqRecord<Bytes>)>, FastqParseError>
    where
        'outer: 'inner1,
        'outer: 'inner2,
    {
        match (
            self.reader1.read_head(container1)?,
            self.reader2.read_head(container2)?,
        ) {
            (Some(()), None) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader2.source.label(),
                    continueal_label: self.reader1.source.label(),
                    eof_pos: self.reader2.source.get_offset(),
                    continueal_pos: self.reader1.source.get_offset() - 1,
                });
            }
            (None, Some(())) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader1.source.label(),
                    continueal_label: self.reader2.source.label(),
                    eof_pos: self.reader1.source.get_offset(),
                    continueal_pos: self.reader2.source.get_offset() - 1,
                });
            }
            (Some(()), Some(())) => {
                let id1 = unsafe { container1.id().unwrap_unchecked() };
                let id2 = unsafe { container2.id().unwrap_unchecked() };
                if id1 != id2 {
                    return Err(FastqParseError::FastqPairError {
                        read1_label: self.reader1.source.label(),
                        read2_label: self.reader2.source.label(),
                        read1_id: String::from_utf8_lossy(id1).into_owned(),
                        read2_id: String::from_utf8_lossy(id2).into_owned(),
                        read1_pos: self.reader1.source.get_offset() - 1,
                        read2_pos: self.reader2.source.get_offset() - 1,
                    });
                }
                self.reader1.read_tail(container1)?;
                self.reader2.read_tail(container2)?;
                return Ok(Some((
                    self.reader1.build(container1),
                    self.reader2.build(container2),
                )));
            }
            (None, None) => {
                return Ok(None);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;

    use tempfile::tempdir;

    use super::*;

    fn to_string(bytes: &[u8]) -> String {
        String::from_utf8(bytes.to_vec()).unwrap()
    }

    #[test]
    fn test_single_end_fastq_chunk_reader() -> Result<()> {
        let tmp = tempdir()?;
        let mut path = PathBuf::from(tmp.path());
        path.push("test.fq");

        let mut file = File::create(&path)?;
        writeln!(
            file,
            "@SEQ_ID1\nGATTA\n+\n!!!!!\n@SEQ_ID2\nACGTA\n+\n#####\n"
        )?;

        let mut chunk_reader = FastqBytesChunkReader::new(BytesReader::new(File::open(&path)?));
        let reader = chunk_reader.chunk_reader()?.unwrap();
        let mut container = FastqContainer::new();
        let record1 = reader.read_record(&mut container).unwrap().unwrap();
        assert_eq!(to_string(&record1.id), "SEQ_ID1");
        assert_eq!(record1.desc, None);
        assert_eq!(to_string(&record1.seq), "GATTA");

        let record2 = reader.read_record(&mut container).unwrap().unwrap();
        assert_eq!(to_string(&record2.id), "SEQ_ID2");
        assert_eq!(record2.desc, None);
        assert_eq!(to_string(&record2.seq), "ACGTA");

        assert!(reader.read_record(&mut container).unwrap().is_none());

        Ok(())
    }

    #[test]
    fn test_paired_end_fastq_chunk_reader() -> Result<()> {
        let tmp = tempdir()?;
        let mut path1 = PathBuf::from(tmp.path());
        let mut path2 = PathBuf::from(tmp.path());
        path1.push("test_1.fq");
        path2.push("test_2.fq");

        let mut f1 = File::create(&path1)?;
        let mut f2 = File::create(&path2)?;

        writeln!(f1, "@SEQ_ID1\nGATTA\n+\n!!!!!\n@SEQ_ID2\nACGTA\n+\n#####")?;
        writeln!(f2, "@SEQ_ID1\nTTTAA\n+\n%%%%%\n@SEQ_ID2\nGGCCC\n+\n&&&&&")?;

        let mut paired_reader = FastqBytesChunkPairedReader::new(
            BytesReader::new(File::open(&path1)?),
            BytesReader::new(File::open(&path2)?),
        );

        let chunk = paired_reader.chunk_reader()?.unwrap();
        let mut container1 = FastqContainer::new();
        let mut container2 = FastqContainer::new();

        let (rec1_r1, rec1_r2) = chunk
            .read_record(&mut container1, &mut container2)?
            .unwrap();
        assert_eq!(to_string(&rec1_r1.id), "SEQ_ID1");
        assert_eq!(to_string(&rec1_r1.seq), "GATTA");
        assert_eq!(to_string(&rec1_r2.id), "SEQ_ID1");
        assert_eq!(to_string(&rec1_r2.seq), "TTTAA");

        let (rec2_r1, rec2_r2) = chunk
            .read_record(&mut container1, &mut container2)?
            .unwrap();
        assert_eq!(to_string(&rec2_r1.id), "SEQ_ID2");
        assert_eq!(to_string(&rec2_r1.seq), "ACGTA");
        assert_eq!(to_string(&rec2_r2.id), "SEQ_ID2");
        assert_eq!(to_string(&rec2_r2.seq), "GGCCC");

        assert!(chunk
            .read_record(&mut container1, &mut container2)?
            .is_none());

        Ok(())
    }
}
