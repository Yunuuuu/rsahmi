use memchr::memchr_iter;

use crate::kractor::reads::parser::fasta::FastaRecord;
use crate::kractor::reads::parser::fastq::{FastqContainer, FastqParseError};

#[derive(Debug)]
pub(crate) struct SliceFastqChunkSource<'a> {
    pos: usize,
    newlines: Vec<usize>,
    newlines_pos: usize,
    chunk: &'a [u8],
    offset: usize,
    label: Option<&'static str>,
}

impl<'a> SliceFastqChunkSource<'a> {
    fn new(chunk: &'a [u8], newlines: Vec<usize>) -> Self {
        Self {
            pos: 0,
            newlines,
            newlines_pos: 0,
            chunk,
            offset: 0,
            label: None,
        }
    }
    fn read_line(&mut self) -> Option<&'a [u8]> {
        if self.newlines_pos >= self.newlines.len() {
            return None;
        }
        let end = unsafe { *self.newlines.get_unchecked(self.newlines_pos) };
        let out = unsafe { self.chunk.get_unchecked(self.pos .. end) };
        self.newlines_pos += 1;
        self.offset += 1;
        self.pos = end + 1;
        Some(out)
    }

    pub fn set_label(&mut self, label: &'static str) {
        self.label = Some(label);
    }

    #[allow(dead_code)]
    pub fn remove_label(&mut self) {
        self.label = None;
    }

    // when split files into multiple chunks, we need adjust position manually.
    #[allow(dead_code)]
    pub fn set_offset(&mut self, offset: usize) {
        self.offset = offset;
    }
}

#[derive(Debug)]
pub(crate) struct MmapFastqReader<'a> {
    source: SliceFastqChunkSource<'a>,
    container: FastqContainer<'a>,
}

impl<'a> MmapFastqReader<'a> {
    pub fn new(source: SliceFastqChunkSource<'a>) -> Self {
        Self {
            source,
            container: FastqContainer::default(),
        }
    }

    pub fn read_record(&mut self) -> Result<Option<FastaRecord<&'a [u8]>>, FastqParseError> {
        // Try reading the 1st line (head). If EOF, return None.
        if self.read_head()?.is_none() {
            return Ok(None);
        };
        self.read_tail()?;
        // SAFETY: We've guaranteed above that 4 lines have been read.
        // This implies the parser has advanced to the Qual state.
        // So `self.build()` will return Some.
        Ok(Some(self.build()))
    }

    fn read_head(&mut self) -> Result<Option<()>, FastqParseError> {
        // Try reading the 1st line (head). If EOF, return None.
        loop {
            match self.source.read_line() {
                Some(line) => {
                    // Skip empty or all-whitespace lines
                    if line.is_empty() || line.iter().all(|b| b.is_ascii_whitespace()) {
                        continue;
                    }

                    self.container
                        .parse_head(line, self.source.label, self.source.offset - 1)?;
                    return Ok(Some(()));
                }
                None => return Ok(None), // EOF
            }
        }
    }

    fn read_tail(&mut self) -> Result<(), FastqParseError> {
        // 2nd line (sequence) must exist. Otherwise, incomplete record.
        if let Some(line) = self.source.read_line() {
            self.container
                .parse_seq(line, self.source.label, self.source.offset - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label,
                record: self.container.to_string(),
                pos: self.source.offset,
            });
        }

        // 3rd line (separator). Must exist.
        if let Some(line) = self.source.read_line() {
            self.container
                .parse_sep(line, self.source.label, self.source.offset - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label,
                record: self.container.to_string(),
                pos: self.source.offset,
            });
        }

        // 4th line (quality). Must exist.
        if let Some(line) = self.source.read_line() {
            self.container
                .parse_qual(line, self.source.label, self.source.offset - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label,
                record: self.container.to_string(),
                pos: self.source.offset,
            });
        }
        Ok(())
    }

    fn build(&mut self) -> FastaRecord<&'a [u8]> {
        let out = unsafe {
            FastaRecord::new(
                self.container.id().unwrap_unchecked(),
                self.container.desc().unwrap_unchecked(),
                self.container.seq().unwrap_unchecked(),
            )
        };
        self.container.reset();
        out
    }
}

pub(crate) struct FastqPairedReader<'a, 'b> {
    reader1: MmapFastqReader<'a>,
    reader2: MmapFastqReader<'b>,
}

impl<'a, 'b> FastqPairedReader<'a, 'b> {
    pub fn new(reader1: MmapFastqReader<'a>, reader2: MmapFastqReader<'b>) -> Self {
        Self { reader1, reader2 }
    }

    pub fn read_record(
        &mut self,
    ) -> Result<Option<(FastaRecord<&'a [u8]>, FastaRecord<&'b [u8]>)>, FastqParseError> {
        match (self.reader1.read_head()?, self.reader2.read_head()?) {
            (Some(()), None) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader2.source.label,
                    continueal_label: self.reader1.source.label,
                    eof_pos: self.reader2.source.offset,
                    continueal_pos: self.reader1.source.offset - 1,
                });
            }
            (None, Some(())) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader1.source.label,
                    continueal_label: self.reader2.source.label,
                    eof_pos: self.reader1.source.offset,
                    continueal_pos: self.reader2.source.offset - 1,
                });
            }
            (Some(()), Some(())) => {
                let id1 = unsafe { self.reader1.container.id().unwrap_unchecked() };
                let id2 = unsafe { self.reader2.container.id().unwrap_unchecked() };
                if id1 != id2 {
                    return Err(FastqParseError::FastqPairError {
                        read1_label: self.reader1.source.label,
                        read2_label: self.reader2.source.label,
                        read1_id: String::from_utf8_lossy(id1).into_owned(),
                        read2_id: String::from_utf8_lossy(id2).into_owned(),
                        read1_pos: self.reader1.source.offset - 1,
                        read2_pos: self.reader2.source.offset - 1,
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

pub(crate) struct SliceChunkReader<'a> {
    label: Option<&'static str>,
    pos: usize,
    line_offset: usize,
    slice: &'a [u8],
    chunk_size: usize,
}

impl<'a> Iterator for SliceChunkReader<'a> {
    type Item = MmapFastqReader<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        self.chunk_reader()
    }
}

impl<'a> SliceChunkReader<'a> {
    #[allow(dead_code)]
    pub fn new(slice: &'a [u8]) -> Self {
        Self::with_capacity(8 * 1024, slice)
    }

    pub fn with_capacity(capacity: usize, slice: &'a [u8]) -> Self {
        Self {
            label: None,
            pos: 0,
            line_offset: 0,
            slice,
            chunk_size: capacity,
        }
    }

    pub fn set_label(&mut self, label: &'static str) {
        self.label = Some(label);
    }

    #[allow(dead_code)]
    pub fn remove_label(&mut self) {
        self.label = None;
    }

    pub fn chunk_reader(&mut self) -> Option<MmapFastqReader<'a>> {
        // If we've reached or passed the end of the input buffer, stop reading
        if self.pos >= self.slice.len() {
            return None;
        }
        let (mut chunk, mut newlines) = slice_chunk(self.slice, self.pos, self.chunk_size);

        // If the chunk contains fewer than 4 lines, the FASTQ record is incomplete
        // so we expand the buffer and try again with a larger chunk
        // we also make sure the file contain other data
        if newlines.len() < 4 && (self.pos + self.chunk_size) < self.slice.len() {
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
        if endpoint == chunk.len() {
            chunk = &chunk[.. endpoint]; // we have no endpoint '\n'
        } else {
            chunk = &chunk[..= endpoint];
        }

        self.pos += chunk.len();
        let line_offset = self.line_offset;
        self.line_offset += newlines.len(); // Increment the number of chunk lines
        let mut source = SliceFastqChunkSource::new(chunk, newlines);
        source.set_offset(line_offset);
        if let Some(label) = self.label {
            source.set_label(label);
        }
        return Some(MmapFastqReader::new(source));
    }
}

pub(crate) struct SliceChunkPairedReader<'a, 'b> {
    label1: Option<&'static str>,
    label2: Option<&'static str>,
    pos1: usize,
    pos2: usize,
    line_offset1: usize,
    line_offset2: usize,
    slice1: &'a [u8],
    slice2: &'b [u8],
    chunk_size: usize,
}

impl<'a, 'b> Iterator for SliceChunkPairedReader<'a, 'b> {
    type Item = std::result::Result<FastqPairedReader<'a, 'b>, FastqParseError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.chunk_reader().transpose()
    }
}

impl<'a, 'b> SliceChunkPairedReader<'a, 'b> {
    #[allow(dead_code)]
    pub fn new(slice1: &'a [u8], slice2: &'b [u8]) -> Self {
        Self::with_capacity(8 * 1024, slice1, slice2)
    }

    pub fn with_capacity(capacity: usize, slice1: &'a [u8], slice2: &'b [u8]) -> Self {
        Self {
            label1: None,
            label2: None,
            pos1: 0,
            pos2: 0,
            line_offset1: 0,
            line_offset2: 0,
            slice1,
            slice2,
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

    pub fn chunk_reader(&mut self) -> Result<Option<FastqPairedReader<'a, 'b>>, FastqParseError> {
        // Handle EOF logic for paired files
        let eof1 = self.pos1 >= self.slice1.len();
        let eof2 = self.pos2 >= self.slice2.len();

        match (eof1, eof2) {
            (true, true) => return Ok(None), // Both inputs exhausted: normal EOF

            // One file is fully read, but the other isn't
            (true, false) => {
                // Allow trailing whitespace in second slice
                if self.slice2[self.pos2 ..]
                    .iter()
                    .all(|b| b.is_ascii_whitespace())
                {
                    self.pos2 = self.slice2.len();
                    return Ok(None);
                }
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.label1,
                    continueal_label: self.label2,
                    eof_pos: self.line_offset1,
                    continueal_pos: self.line_offset2,
                });
            }
            (false, true) => {
                if self.slice1[self.pos1 ..]
                    .iter()
                    .all(|b| b.is_ascii_whitespace())
                {
                    self.pos1 = self.slice1.len();
                    return Ok(None);
                }
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.label2,
                    continueal_label: self.label1,
                    eof_pos: self.line_offset2,
                    continueal_pos: self.line_offset1,
                });
            }
            (false, false) => {} // Continue reading normally
        }
        let end1 = self.pos1 + self.chunk_size;
        let end2 = self.pos2 + self.chunk_size;

        // Slice each input at current position to retrieve a chunk and its line breaks
        let (mut chunk1, mut newlines1) = slice_chunk(self.slice1, self.pos1, self.chunk_size);
        let (mut chunk2, mut newlines2) = slice_chunk(self.slice2, self.pos2, self.chunk_size);

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        let count = usize::min(newlines1.len(), newlines2.len());

        // If we didn't capture a complete record
        if count < 4 {
            // more data is available,
            // double chunk size and try again (adaptive chunking)
            if end1 < self.slice1.len() && end2 < self.slice2.len() {
                self.chunk_size *= 2;
                return self.chunk_reader();
            } else {
                // If one side has incomplete records at EOF, raise error
                let (label, pos, record) = if end1 >= self.slice1.len() {
                    (
                        self.label1,
                        self.line_offset1,
                        String::from_utf8_lossy(chunk1).to_string(),
                    )
                } else {
                    (
                        self.label2,
                        self.line_offset2,
                        String::from_utf8_lossy(chunk2).to_string(),
                    )
                };
                return Err(FastqParseError::IncompleteRecord { label, record, pos });
            }
        }

        // Truncate to a multiple of 4 lines to ensure whole records only
        let nkeep = count - (count % 4);
        newlines1.truncate(nkeep);
        newlines2.truncate(nkeep);

        // SAFETY: newlines won't be empty
        let endpoint1 = unsafe { *newlines1.last().unwrap_unchecked() };
        // no trailing newline
        if endpoint1 < chunk1.len() {
            chunk1 = &chunk1[..= endpoint1]; // include newline
        }
        let endpoint2 = unsafe { *newlines2.last().unwrap_unchecked() };
        // no trailing newline
        if endpoint2 < chunk2.len() {
            chunk2 = &chunk2[..= endpoint2];
        }

        // Update position offsets for next chunk
        self.pos1 += chunk1.len();
        self.pos2 += chunk2.len();

        // Line offsets are in terms of logical FASTQ lines
        let line_offset1 = self.line_offset1;
        let line_offset2 = self.line_offset2;

        self.line_offset1 += newlines1.len();
        self.line_offset2 += newlines2.len();

        // Wrap chunks into source readers with correct offset and label
        let mut source1 = SliceFastqChunkSource::new(chunk1, newlines1);
        source1.set_offset(line_offset1);
        if let Some(label1) = self.label1 {
            source1.set_label(label1);
        }
        let mut source2 = SliceFastqChunkSource::new(chunk2, newlines2);
        source2.set_offset(line_offset2);
        if let Some(label2) = self.label2 {
            source2.set_label(label2);
        }
        // Return paired reader
        Ok(Some(FastqPairedReader::new(
            MmapFastqReader::new(source1),
            MmapFastqReader::new(source2),
        )))
    }
}

fn slice_chunk<'a>(slice: &'a [u8], pos: usize, chunk_size: usize) -> (&'a [u8], Vec<usize>) {
    let chunk;
    let mut newlines: Vec<usize>;
    let end = pos + chunk_size;
    // take consideration of the last chunk
    if end >= slice.len() {
        // SAFETY: We're within bounds because pos < bytes.len()
        chunk = *unsafe { &slice.get_unchecked(pos ..) };
        // Find all newline character positions in the chunk
        newlines = memchr_iter(b'\n', chunk).collect();

        // If there are no newlines, treat the whole chunk as a single (newline-less) line
        if newlines.is_empty() {
            newlines.push(chunk.len());
        } else if (chunk.len() - 1)
            // Ensure the last line is included even if it doesn't end with `\n`
            // SAFETY: newlines is non-empty, so unwrap_unchecked is safe
                != *unsafe { newlines.get_unchecked(newlines.len() - 1) }
        {
            newlines.push(chunk.len());
        }
    } else {
        // SAFETY: end < bytes.len() is guaranteed here
        chunk = *unsafe { &slice.get_unchecked(pos .. end) };

        // Find newline offsets within this chunk
        newlines = memchr_iter(b'\n', chunk).collect();
    }
    (chunk, newlines)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTQ_SINGLE: &str =
        "@SEQ_ID1\nGATTACA\n+\n!!!!!!!\n@SEQ_ID2\nCTGAAGT\n+\n???????\n@SEQ_ID3\nTTTGGCA\n+\n%%%%%%%\n@SEQ_ID4\nAACCGGT\n+\n(((((((\n";
    const FASTQ_PAIR_1: &str =
        "@SEQ_ID1\nGATTACA\n+\n!!!!!!!\n@SEQ_ID2\nCTGAAGT\n+\n???????\n@SEQ_ID3\nTTTGGCA\n+\n%%%%%%%\n@SEQ_ID4\nAACCGGT\n+\n(((((((\n";

    const FASTQ_PAIR_2: &str =
        "@SEQ_ID1\nTTGGAAC\n+\n!!!!!!!\n@SEQ_ID2\nACTTCAG\n+\n???????\n@SEQ_ID3\nTGCCAAA\n+\n%%%%%%%\n@SEQ_ID4\nACCGGTT\n+\n(((((((\n";

    #[test]
    fn test_slice_chunk_reader_single() {
        let bytes = FASTQ_SINGLE.as_bytes();
        let reader = SliceChunkReader::with_capacity(32, bytes);

        let mut total_records = 0;
        for mut chunk in reader {
            while let Some(record) = chunk.read_record().unwrap() {
                assert!(record.seq.len() > 0);
                total_records += 1;
            }
        }
        assert_eq!(total_records, 4);
    }

    #[test]
    fn test_slice_chunk_reader_paired() {
        let fq1 = FASTQ_PAIR_1.as_bytes();
        let fq2 = FASTQ_PAIR_2.as_bytes();
        let reader = SliceChunkPairedReader::with_capacity(32, fq1, fq2);

        let mut total_records = 0;
        for pair in reader {
            let mut pair = pair.unwrap();
            while let Some((r1, r2)) = pair.read_record().unwrap() {
                assert_eq!(r1.id[5 ..], r2.id[5 ..]); // match SEQ_IDn
                total_records += 1;
            }
        }
        assert_eq!(total_records, 4);
    }
}
