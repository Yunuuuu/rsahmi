use crate::parser::fastq::{container::FastqContainer, FastqParseError, FastqRecord};
use crate::reader::slice::{SliceLineReader, SliceProgressBarReader};

pub(crate) struct FastqSliceChunkReader<'a> {
    line_offset: usize,
    chunk_size: usize,
    reader: SliceProgressBarReader<'a>,
}

impl<'a> FastqSliceChunkReader<'a> {
    #[allow(dead_code)]
    pub(crate) fn new(reader: SliceProgressBarReader<'a>) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub(crate) fn with_capacity(capacity: usize, reader: SliceProgressBarReader<'a>) -> Self {
        Self {
            line_offset: 0,
            chunk_size: capacity,
            reader,
        }
    }

    pub(crate) fn chunk_reader(
        &mut self,
    ) -> Option<FastqSliceReader<'a, std::vec::IntoIter<usize>>> {
        // If we've reached or passed the end of the input buffer, stop reading
        let slice = match self.reader.take_slice(self.chunk_size) {
            Some(chunk) => chunk,
            None => return None,
        };
        let mut newlines = slice.newlines();
        let mut chunk = slice.as_slice();

        // If the chunk contains fewer than 4 lines, the FASTQ record is incomplete
        // so we expand the buffer and try again with a larger chunk
        // we also make sure the file contain other data
        if newlines.len() < 4 && !slice.eof() {
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
        self.reader.advance(chunk.len());

        let line_offset = self.line_offset;
        self.line_offset += newlines.len(); // Increment the number of chunk lines
        let mut source = FastqSliceSource::new(chunk, newlines.into_iter());
        source.set_offset(line_offset);
        if let Some(label) = self.reader.label() {
            source.set_label(label);
        }
        return Some(FastqSliceReader::new(source));
    }
}

pub(crate) struct FastqSliceChunkPairedReader<'a, 'b> {
    line_offset1: usize,
    line_offset2: usize,
    chunk_size: usize,
    reader1: SliceProgressBarReader<'a>,
    reader2: SliceProgressBarReader<'b>,
}

impl<'a, 'b> FastqSliceChunkPairedReader<'a, 'b> {
    #[allow(dead_code)]
    pub(crate) fn new(
        reader1: SliceProgressBarReader<'a>,
        reader2: SliceProgressBarReader<'b>,
    ) -> Self {
        Self::with_capacity(8 * 1024, reader1, reader2)
    }

    pub(crate) fn with_capacity(
        capacity: usize,
        reader1: SliceProgressBarReader<'a>,
        reader2: SliceProgressBarReader<'b>,
    ) -> Self {
        Self {
            line_offset1: 0,
            line_offset2: 0,
            chunk_size: capacity,
            reader1: reader1,
            reader2: reader2,
        }
    }

    pub(crate) fn chunk_reader(
        &mut self,
    ) -> Result<Option<FastqSlicePairedReader<'a, 'b, std::vec::IntoIter<usize>>>, FastqParseError>
    {
        // Handle EOF logic for paired files
        let eof1 = self.reader1.eof();
        let eof2 = self.reader2.eof();

        match (eof1, eof2) {
            (true, true) => return Ok(None), // Both inputs exhausted: normal EOF

            // One file is fully read, but the other isn't
            (true, false) => {
                // Allow trailing whitespace in second slice
                if self
                    .reader2
                    .as_slice()
                    .iter()
                    .all(|b| b.is_ascii_whitespace())
                {
                    self.reader2.finish();
                    return Ok(None);
                }
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader2.label(),
                    continueal_label: self.reader1.label(),
                    eof_pos: self.line_offset1,
                    continueal_pos: self.line_offset2,
                });
            }
            (false, true) => {
                if self
                    .reader1
                    .as_slice()
                    .iter()
                    .all(|b| b.is_ascii_whitespace())
                {
                    self.reader1.finish();
                    return Ok(None);
                }
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader2.label(),
                    continueal_label: self.reader1.label(),
                    eof_pos: self.line_offset2,
                    continueal_pos: self.line_offset1,
                });
            }
            (false, false) => {} // Continue reading normally
        }
        let slice1 = unsafe { self.reader1.take_slice(self.chunk_size).unwrap_unchecked() };
        let slice2 = unsafe { self.reader2.take_slice(self.chunk_size).unwrap_unchecked() };

        // Slice each input at current position to retrieve a chunk and its line breaks
        let mut newlines1 = slice1.newlines();
        let mut chunk1 = slice1.as_slice();
        let mut newlines2 = slice2.newlines();
        let mut chunk2 = slice2.as_slice();

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        let count = usize::min(newlines1.len(), newlines2.len());

        // If we didn't capture a complete record
        if count < 4 {
            // more data is available,
            // double chunk size and try again (adaptive chunking)
            if !slice1.eof() && !slice2.eof() {
                self.chunk_size *= 2;
                return self.chunk_reader();
            } else {
                // If one side has incomplete records at EOF, raise error
                let (label, pos, record) = if slice1.eof() {
                    (
                        self.reader1.label(),
                        self.line_offset1,
                        String::from_utf8_lossy(chunk1).to_string(),
                    )
                } else {
                    (
                        self.reader2.label(),
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
        self.reader1.advance(chunk1.len());
        self.reader2.advance(chunk2.len());

        // Line offsets are in terms of logical FASTQ lines
        let line_offset1 = self.line_offset1;
        let line_offset2 = self.line_offset2;

        self.line_offset1 += newlines1.len();
        self.line_offset2 += newlines2.len();

        // Wrap chunks into source readers with correct offset and label
        let mut source1 = FastqSliceSource::new(chunk1, newlines1.into_iter());

        source1.set_offset(line_offset1);
        if let Some(label1) = self.reader1.label() {
            source1.set_label(label1);
        }
        let mut source2 = FastqSliceSource::new(chunk2, newlines2.into_iter());
        source2.set_offset(line_offset2);
        if let Some(label2) = self.reader2.label() {
            source2.set_label(label2);
        }

        // Return paired reader
        Ok(Some(FastqSlicePairedReader::new(
            FastqSliceReader::new(source1),
            FastqSliceReader::new(source2),
        )))
    }
}

#[derive(Debug)]
pub(crate) struct FastqSliceSource<'a, I: Iterator<Item = usize>> {
    offset: usize,
    reader: SliceLineReader<'a, I>,
}

impl<'a, I: Iterator<Item = usize>> FastqSliceSource<'a, I> {
    fn new(slice: &'a [u8], breaks: I) -> Self {
        Self {
            offset: 0,
            reader: SliceLineReader::with_breaks(breaks, slice),
        }
    }

    fn read_line(&mut self) -> Option<&'a [u8]> {
        let out = self.reader.read_line()?;
        self.offset += 1;
        Some(out)
    }

    fn set_label(&mut self, label: &'static str) {
        self.reader.set_label(label)
    }

    fn label(&self) -> Option<&'static str> {
        self.reader.label()
    }

    // when split files into multiple chunks, we need adjust position manually.
    #[allow(dead_code)]
    pub fn set_offset(&mut self, offset: usize) {
        self.offset = offset;
    }
}

#[derive(Debug)]
pub(crate) struct FastqSliceReader<'a, I: Iterator<Item = usize>> {
    source: FastqSliceSource<'a, I>,
    container: FastqContainer<'a>,
}

impl<'a, I: Iterator<Item = usize>> FastqSliceReader<'a, I> {
    pub(crate) fn new(source: FastqSliceSource<'a, I>) -> Self {
        Self {
            source,
            container: FastqContainer::default(),
        }
    }

    pub(crate) fn read_record(&mut self) -> Result<Option<FastqRecord<&'a [u8]>>, FastqParseError> {
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
                        .parse_head(line, self.source.label(), self.source.offset - 1)?;
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
                .parse_seq(line, self.source.label(), self.source.offset - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label(),
                record: self.container.to_string(),
                pos: self.source.offset,
            });
        }

        // 3rd line (separator). Must exist.
        if let Some(line) = self.source.read_line() {
            self.container
                .parse_sep(line, self.source.label(), self.source.offset - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label(),
                record: self.container.to_string(),
                pos: self.source.offset,
            });
        }

        // 4th line (quality). Must exist.
        if let Some(line) = self.source.read_line() {
            self.container
                .parse_qual(line, self.source.label(), self.source.offset - 1)?;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.source.label(),
                record: self.container.to_string(),
                pos: self.source.offset,
            });
        }
        Ok(())
    }

    fn build(&mut self) -> FastqRecord<&'a [u8]> {
        let out = unsafe {
            FastqRecord::new(
                self.container.id().unwrap_unchecked(),
                self.container.desc().unwrap_unchecked(),
                self.container.seq().unwrap_unchecked(),
                self.container.sep().unwrap_unchecked(),
                self.container.qual().unwrap_unchecked(),
            )
        };
        self.container.reset();
        out
    }
}

pub(crate) struct FastqSlicePairedReader<'a, 'b, I: Iterator<Item = usize>> {
    reader1: FastqSliceReader<'a, I>,
    reader2: FastqSliceReader<'b, I>,
}

impl<'a, 'b, I: Iterator<Item = usize>> FastqSlicePairedReader<'a, 'b, I> {
    pub(crate) fn new(reader1: FastqSliceReader<'a, I>, reader2: FastqSliceReader<'b, I>) -> Self {
        Self { reader1, reader2 }
    }

    pub(crate) fn read_record(
        &mut self,
    ) -> Result<Option<(FastqRecord<&'a [u8]>, FastqRecord<&'b [u8]>)>, FastqParseError> {
        match (self.reader1.read_head()?, self.reader2.read_head()?) {
            (Some(()), None) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader2.source.label(),
                    continueal_label: self.reader1.source.label(),
                    eof_pos: self.reader2.source.offset,
                    continueal_pos: self.reader1.source.offset - 1,
                });
            }
            (None, Some(())) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader1.source.label(),
                    continueal_label: self.reader2.source.label(),
                    eof_pos: self.reader1.source.offset,
                    continueal_pos: self.reader2.source.offset - 1,
                });
            }
            (Some(()), Some(())) => {
                let id1 = unsafe { self.reader1.container.id().unwrap_unchecked() };
                let id2 = unsafe { self.reader2.container.id().unwrap_unchecked() };
                if id1 != id2 {
                    return Err(FastqParseError::FastqPairError {
                        read1_label: self.reader1.source.label(),
                        read2_label: self.reader2.source.label(),
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

        let mut reader =
            FastqSliceChunkReader::with_capacity(32, SliceProgressBarReader::new(bytes));

        let mut total_records = 0;
        while let Some(mut chunk) = reader.chunk_reader() {
            while let Some(record) = chunk.read_record().unwrap() {
                assert!(record.seq.len() > 0);
                total_records += 1;
            }
        }
        assert_eq!(total_records, 4);
    }

    #[test]
    fn test_slice_chunk_reader_paired() {
        let fq1 = SliceProgressBarReader::new(FASTQ_PAIR_1.as_bytes());
        let fq2 = SliceProgressBarReader::new(FASTQ_PAIR_2.as_bytes());
        let mut reader = FastqSliceChunkPairedReader::with_capacity(32, fq1, fq2);

        let mut total_records = 0;
        while let Some(mut pair_chunk) = reader.chunk_reader().unwrap() {
            while let Some((r1, r2)) = pair_chunk.read_record().unwrap() {
                assert_eq!(r1.id[5 ..], r2.id[5 ..]); // match SEQ_IDn
                total_records += 1;
            }
        }
        assert_eq!(total_records, 4);
    }
}
