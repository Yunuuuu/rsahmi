use std::cell::Cell;
use std::cell::RefCell;
use std::io::Read;

use anyhow::{anyhow, Result};
use bytes::{Bytes, BytesMut};
use memchr::memchr_iter;

use crate::kractor::reads::parser::fasta::FastaRecord;
use crate::kractor::reads::parser::fastq::{FastqContainer, FastqParseError};

#[derive(Debug)]
pub(crate) struct BytesFastqChunkSource {
    pos: Cell<usize>,
    newlines: Vec<usize>,
    newlines_pos: Cell<usize>,
    // we need hold long-lived references to the chunk, so we use
    // normal reference, and wrap pos and newlines_pos into Cell
    // so we can change them
    chunk: Bytes,
    offset: Cell<usize>,
    label: Option<&'static str>,
}

impl BytesFastqChunkSource {
    fn new(chunk: Bytes, newlines: Vec<usize>) -> Self {
        Self {
            pos: Cell::new(0),
            newlines,
            newlines_pos: Cell::new(0),
            chunk: chunk,
            offset: Cell::new(0),
            label: None,
        }
    }

    fn read_line(&self) -> Option<&[u8]> {
        if self.newlines_pos.get() >= self.newlines.len() {
            return None;
        }
        let end = unsafe { *self.newlines.get_unchecked(self.newlines_pos.get()) };
        let out = unsafe { self.chunk.get_unchecked(self.pos.get() .. end) };

        self.newlines_pos.set(self.newlines_pos.get() + 1);
        self.offset.set(self.offset.get() + 1);
        self.pos.set(end + 1);
        Some(out)
    }

    fn offset(&self) -> usize {
        self.offset.get()
    }

    fn label(&self) -> Option<&'static str> {
        self.label
    }

    fn set_label(&mut self, label: &'static str) {
        self.label = Some(label);
    }

    #[allow(dead_code)]
    fn remove_label(&mut self) {
        self.label = None;
    }

    // when split files into multiple chunks, we need adjust position manually.
    fn set_offset(&mut self, offset: usize) {
        self.offset = Cell::new(offset);
    }
}

#[derive(Debug)]
pub(crate) struct BytesFastqReader<'a> {
    source: BytesFastqChunkSource,
    container: RefCell<FastqContainer<'a>>,
}

impl<'s, 'a> BytesFastqReader<'a>
where
    's: 'a,
{
    pub fn new(source: BytesFastqChunkSource) -> Self {
        Self {
            source,
            container: RefCell::new(FastqContainer::<'a>::default()),
        }
    }

    fn read_record(&'s self) -> Result<Option<FastaRecord<Bytes>>, FastqParseError> {
        // Try reading the 1st line (head). If EOF, return None.
        if self.read_head()?.is_none() {
            return Ok(None);
        };
        self.read_tail()?;
        Ok(Some(self.build()))
    }

    fn read_head(&'s self) -> Result<Option<()>, FastqParseError> {
        // Try reading the 1st line (head). If EOF, return None.
        loop {
            match self.source.read_line() {
                Some(line) => {
                    // Skip empty or all-whitespace lines
                    if line.is_empty() || line.iter().all(|b| b.is_ascii_whitespace()) {
                        continue;
                    }
                    self.container.borrow_mut().parse_head(
                        line,
                        self.source.label(),
                        self.source.offset() - 1,
                    )?;
                    return Ok(Some(()));
                }
                None => return Ok(None), // EOF
            }
        }
    }

    fn read_tail(&'s self) -> Result<(), FastqParseError> {
        // 2nd line (sequence) must exist. Otherwise, incomplete record.
        if let Some(line) = self.source.read_line() {
            self.container.borrow_mut().parse_seq(
                line,
                self.source.label,
                self.source.offset() - 1,
            )?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label,
                record: self.container.borrow().to_string(),
                pos: self.source.offset(),
            });
        }

        // 3rd line (separator). Must exist.
        if let Some(line) = self.source.read_line() {
            self.container.borrow_mut().parse_sep(
                line,
                self.source.label,
                self.source.offset() - 1,
            )?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label,
                record: self.container.borrow().to_string(),
                pos: self.source.offset(),
            });
        }

        // 4th line (quality). Must exist.
        if let Some(line) = self.source.read_line() {
            self.container.borrow_mut().parse_qual(
                line,
                self.source.label,
                self.source.offset() - 1,
            )?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label,
                record: self.container.borrow().to_string(),
                pos: self.source.offset(),
            });
        }
        Ok(())
    }

    fn build(&self) -> FastaRecord<Bytes> {
        let out = unsafe {
            FastaRecord::new(
                self.source
                    .chunk
                    .slice_ref(self.container.borrow().id().unwrap_unchecked()),
                self.container
                    .borrow()
                    .desc()
                    .unwrap_unchecked()
                    .map(|slice| self.source.chunk.slice_ref(slice)),
                self.source
                    .chunk
                    .slice_ref(self.container.borrow().seq().unwrap_unchecked()),
            )
        };
        self.container.borrow_mut().reset();
        out
    }
}

pub(crate) struct BytesFastqPairedChunkReader<'a, 'b> {
    reader1: BytesFastqReader<'a>,
    reader2: BytesFastqReader<'b>,
}

impl<'s, 'a, 'b> BytesFastqPairedChunkReader<'a, 'b>
where
    's: 'a,
    's: 'b,
{
    pub fn new(reader1: BytesFastqReader<'a>, reader2: BytesFastqReader<'b>) -> Self {
        Self { reader1, reader2 }
    }

    pub fn read_record(
        &'s self,
    ) -> Result<Option<(FastaRecord<Bytes>, FastaRecord<Bytes>)>, FastqParseError> {
        match (self.reader1.read_head()?, self.reader2.read_head()?) {
            (Some(()), None) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader2.source.label,
                    continueal_label: self.reader1.source.label,
                    eof_pos: self.reader2.source.offset(),
                    continueal_pos: self.reader1.source.offset() - 1,
                });
            }
            (None, Some(())) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader1.source.label,
                    continueal_label: self.reader2.source.label,
                    eof_pos: self.reader1.source.offset(),
                    continueal_pos: self.reader2.source.offset() - 1,
                });
            }
            (Some(()), Some(())) => {
                let id1 = unsafe { self.reader1.container.borrow().id().unwrap_unchecked() };
                let id2 = unsafe { self.reader2.container.borrow().id().unwrap_unchecked() };
                if id1 != id2 {
                    return Err(FastqParseError::FastqPairError {
                        read1_label: self.reader1.source.label,
                        read2_label: self.reader2.source.label,
                        read1_id: String::from_utf8_lossy(id1).into_owned(),
                        read2_id: String::from_utf8_lossy(id2).into_owned(),
                        read1_pos: self.reader1.source.offset() - 1,
                        read2_pos: self.reader2.source.offset() - 1,
                    });
                }
                self.reader1.read_tail()?;
                self.reader2.read_tail()?;
                return Ok(Some((self.reader1.build(), self.reader2.build())));
            }
            (None, None) => {
                return Ok(None);
            }
        }
    }
}

pub(crate) struct BytesChunkReader<R>
where
    R: Read,
{
    label: Option<&'static str>,
    line_offset: usize,
    reader: R,
    chunk_size: usize,
    leftover: BytesMut, // stores bytes after the last \n
}

impl<R> BytesChunkReader<R>
where
    R: Read,
{
    #[allow(dead_code)]
    pub fn new(reader: R) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Self {
            label: None,
            line_offset: 0,
            reader,
            chunk_size: capacity,
            leftover: BytesMut::new(),
        }
    }

    pub fn set_label(&mut self, label: &'static str) {
        self.label = Some(label);
    }

    #[allow(dead_code)]
    pub fn remove_label(&mut self) {
        self.label = None;
    }

    pub fn chunk_reader<'a>(&'a mut self) -> Result<Option<BytesFastqReader<'a>>> {
        // read files
        let (mut buf, nbytes) = read(&mut self.reader, &self.leftover, self.chunk_size)?;

        // If we've reached or passed the end of the input buffer, stop reading
        if nbytes == 0 {
            self.leftover = BytesMut::new();
            if buf.is_empty() {
                return Ok(None);
            } else {
                // take consideration of the last chunk
                let mut newlines: Vec<usize> = memchr_iter(b'\n', &buf).collect();
                // If there are no newlines, treat the whole chunk as a single (newline-less) line
                if newlines.is_empty() {
                    newlines.push(buf.len());
                } else if (buf.len() - 1)
                // Ensure the last line is included even if it doesn't end with `\n`
                // SAFETY: newlines is non-empty, so unwrap_unchecked is safe
                != *unsafe { newlines.get_unchecked(newlines.len() - 1) }
                {
                    newlines.push(buf.len());
                }
                let line_offset = self.line_offset;
                self.line_offset += newlines.len(); // Increment the number of chunk lines
                let mut source = BytesFastqChunkSource::new(buf.freeze(), newlines);
                source.set_offset(line_offset);
                if let Some(label) = self.label {
                    source.set_label(label);
                }
                return Ok(Some(BytesFastqReader::new(source)));
            }
        }

        // Find all newline character positions in the chunk
        let mut newlines: Vec<usize> = memchr_iter(b'\n', &buf).collect();

        // If the chunk contains fewer than 4 lines, the FASTQ record is incomplete
        // so we expand the buffer and try again with a larger chunk
        if newlines.len() < 4 {
            self.leftover = buf; // remainder for next round
                                 // we always double the buffer size, if there no enough buffer to hold a whole fastq record
            self.chunk_size *= 2;
            return self.chunk_reader();
        }

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        newlines.truncate(newlines.len() - (newlines.len() % 4));

        // SAFETY: above code will ensure newlines won't be empty
        let endpoint = unsafe { *newlines.last().unwrap_unchecked() };
        let chunk = buf.split_to(endpoint + 1);
        self.leftover = buf;

        let line_offset = self.line_offset;
        self.line_offset += newlines.len(); // Increment the number of chunk lines
        let mut source = BytesFastqChunkSource::new(chunk.freeze(), newlines);
        source.set_offset(line_offset);
        if let Some(label) = self.label {
            source.set_label(label);
        }
        return Ok(Some(BytesFastqReader::new(source)));
    }
}

pub(crate) struct BytesChunkPairedReader<R>
where
    R: Read,
{
    label1: Option<&'static str>,
    label2: Option<&'static str>,
    line_offset1: usize,
    line_offset2: usize,
    reader1: R,
    leftover1: BytesMut,
    reader2: R,
    leftover2: BytesMut,
    chunk_size: usize,
}

impl<R> BytesChunkPairedReader<R>
where
    R: Read,
{
    #[allow(dead_code)]
    pub fn new(reader1: R, reader2: R) -> Self {
        Self::with_capacity(8 * 1024, reader1, reader2)
    }

    pub fn with_capacity(capacity: usize, reader1: R, reader2: R) -> Self {
        Self {
            label1: None,
            label2: None,
            line_offset1: 0,
            line_offset2: 0,
            reader1,
            leftover1: BytesMut::new(),
            reader2,
            leftover2: BytesMut::new(),
            chunk_size: capacity,
        }
    }

    pub fn set_label1(&mut self, label: &'static str) {
        self.label1 = Some(label);
    }

    pub fn set_label2(&mut self, label: &'static str) {
        self.label2 = Some(label);
    }

    #[allow(dead_code)]
    pub fn remove_label1(&mut self) {
        self.label1 = None;
    }

    #[allow(dead_code)]
    pub fn remove_label2(&mut self) {
        self.label2 = None;
    }

    pub fn chunk_reader(&mut self) -> Result<Option<BytesFastqPairedChunkReader>> {
        // --- Read data for read1 and append to existing leftover ---
        let (mut buf1, nbytes1) = read(&mut self.reader1, &self.leftover1, self.chunk_size)?;
        let eof1 = nbytes1 == 0 && buf1.is_empty();

        // --- Read data for read2 and append to existing leftover ---
        let (mut buf2, nbytes2) = read(&mut self.reader2, &self.leftover2, self.chunk_size)?;
        let eof2 = nbytes2 == 0 && buf2.is_empty();

        // --- EOF + desync checks ---
        match (eof1, eof2) {
            (true, true) => return Ok(None),

            // read1 ended, check if read2 has only trailing whitespace
            (true, false) => {
                if buf2.iter().all(|b| b.is_ascii_whitespace()) {
                    // read and discard trailing whitespace from stream
                    for byte_result in self.reader2.by_ref().bytes() {
                        let byte = byte_result?;
                        if !byte.is_ascii_whitespace() {
                            break;
                        }
                    }
                    return Ok(None);
                }
                // Possible true sync failure
                return Err(anyhow!("{}", FastqParseError::OutOfSync {
                    eof_label: self.label1,
                    continueal_label: self.label2,
                    eof_pos: self.line_offset1,
                    continueal_pos: self.line_offset2,
                }));
            }
            (false, true) => {
                if buf1.iter().all(|b| b.is_ascii_whitespace()) {
                    for byte_result in self.reader1.by_ref().bytes() {
                        let byte = byte_result?;
                        if !byte.is_ascii_whitespace() {
                            break;
                        }
                    }
                    return Ok(None);
                }
                return Err(anyhow!("{}", FastqParseError::OutOfSync {
                    eof_label: self.label2,
                    continueal_label: self.label1,
                    eof_pos: self.line_offset2,
                    continueal_pos: self.line_offset1,
                }));
            }
            (false, false) => {}
        }

        // --- Find newlines in both chunks ---
        // read1
        let mut newlines1: Vec<usize> = memchr_iter(b'\n', &buf1).collect();
        if nbytes1 == 0 {
            // Edge case: input ends without newline
            // no data left, we ensure the final line is counted
            if newlines1.is_empty() {
                newlines1.push(buf1.len());
            } else if (buf1.len() - 1)
            // Ensure the last line is included even if it doesn't end with `\n`
            // SAFETY: newlines is non-empty, so unwrap_unchecked is safe
            != *unsafe { newlines1.get_unchecked(newlines1.len() - 1) }
            {
                newlines1.push(buf1.len());
            }
        }

        // read2
        let mut newlines2: Vec<usize> = memchr_iter(b'\n', &buf2).collect();
        if nbytes2 == 0 {
            // Edge case: input ends without newline
            // no data left, we ensure the final line is counted
            if newlines2.is_empty() {
                newlines2.push(buf2.len());
            } else if (buf2.len() - 1)
                // Ensure the last line is included even if it doesn't end with `\n`
                // SAFETY: newlines is non-empty, so unwrap_unchecked is safe
                != *unsafe { newlines2.get_unchecked(newlines2.len() - 1) }
            {
                newlines2.push(buf2.len());
            }
        }

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        // --- Determine how many records we can safely return (must be multiple of 4 lines) ---
        let count = usize::min(newlines1.len(), newlines2.len());
        if count < 4 {
            if nbytes1 > 0 && nbytes2 > 0 {
                // Not enough lines yet, double chunk size for next round
                self.leftover1 = buf1;
                self.leftover2 = buf2;
                self.chunk_size *= 2;
                return self.chunk_reader();
            } else {
                // Partial record at EOF — invalid
                // If one file only contain few lines, it means it has incompleted record
                let (label, pos, record) = if nbytes1 == 0 {
                    (
                        self.label1,
                        self.line_offset1,
                        String::from_utf8_lossy(&buf1).to_string(),
                    )
                } else {
                    (
                        self.label2,
                        self.line_offset2,
                        String::from_utf8_lossy(&buf2).to_string(),
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
        let chunk1;
        if endpoint1 < buf1.len() {
            chunk1 = buf1.split_to(endpoint1 + 1); // include `\n`
            self.leftover1 = buf1;
        } else {
            chunk1 = buf1;
            self.leftover1 = BytesMut::new();
        }

        // --- Split buf2 likewise ---
        let endpoint2 = unsafe { *newlines2.last().unwrap_unchecked() };
        let chunk2;
        if endpoint2 < buf2.len() {
            chunk2 = buf2.split_to(endpoint2 + 1); // include `\n`
            self.leftover2 = buf2;
        } else {
            chunk2 = buf2;
            self.leftover2 = BytesMut::new();
        }

        // --- Update line offset bookkeeping ---
        let line_offset1 = self.line_offset1;
        let line_offset2 = self.line_offset2;
        self.line_offset1 += newlines1.len();
        self.line_offset2 += newlines2.len();

        // --- Wrap into sources with offsets and labels ---
        let mut source1 = BytesFastqChunkSource::new(chunk1.freeze(), newlines1);
        source1.set_offset(line_offset1);
        if let Some(label1) = self.label1 {
            source1.set_label(label1);
        }
        let mut source2 = BytesFastqChunkSource::new(chunk2.freeze(), newlines2);
        source2.set_offset(line_offset2);
        if let Some(label2) = self.label2 {
            source2.set_label(label2);
        }
        Ok(Some(BytesFastqPairedChunkReader::new(
            BytesFastqReader::new(source1),
            BytesFastqReader::new(source2),
        )))
    }
}

fn read<R: Read>(
    read: &mut R,
    leftover: &BytesMut,
    chunk_size: usize,
) -> Result<(BytesMut, usize)> {
    let leftover_len = leftover.len();
    let mut buf = BytesMut::with_capacity(chunk_size + leftover_len);
    buf.extend_from_slice(leftover);

    unsafe { buf.set_len(leftover_len + chunk_size) };
    let nbytes = read.read(&mut buf[leftover_len ..])?;
    unsafe { buf.set_len(leftover_len + nbytes) };
    Ok((buf, nbytes))
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

        let mut chunk_reader = BytesChunkReader::new(File::open(&path)?);
        let reader = chunk_reader.chunk_reader()?.unwrap();

        let record1 = reader.read_record().unwrap().unwrap();
        assert_eq!(to_string(&record1.id), "SEQ_ID1");
        assert_eq!(record1.desc, None);
        assert_eq!(to_string(&record1.seq), "GATTA");

        let record2 = reader.read_record().unwrap().unwrap();
        assert_eq!(to_string(&record2.id), "SEQ_ID2");
        assert_eq!(record2.desc, None);
        assert_eq!(to_string(&record2.seq), "ACGTA");

        assert!(reader.read_record().unwrap().is_none());

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

        let mut paired_reader =
            BytesChunkPairedReader::new(File::open(&path1)?, File::open(&path2)?);

        let chunk = paired_reader.chunk_reader()?.unwrap();

        let (rec1_r1, rec1_r2) = chunk.read_record()?.unwrap();
        assert_eq!(to_string(&rec1_r1.id), "SEQ_ID1");
        assert_eq!(to_string(&rec1_r1.seq), "GATTA");
        assert_eq!(to_string(&rec1_r2.id), "SEQ_ID1");
        assert_eq!(to_string(&rec1_r2.seq), "TTTAA");

        let (rec2_r1, rec2_r2) = chunk.read_record()?.unwrap();
        assert_eq!(to_string(&rec2_r1.id), "SEQ_ID2");
        assert_eq!(to_string(&rec2_r1.seq), "ACGTA");
        assert_eq!(to_string(&rec2_r2.id), "SEQ_ID2");
        assert_eq!(to_string(&rec2_r2.seq), "GGCCC");

        assert!(chunk.read_record()?.is_none());

        Ok(())
    }
}
