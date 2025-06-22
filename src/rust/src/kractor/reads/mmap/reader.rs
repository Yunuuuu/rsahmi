use memchr::memchr_iter;

use crate::parser::fasta::FastaRecord;
use crate::parser::fastq::{
    FastqPairedReader, FastqParseError, FastqReader, FastqSource,
};

pub struct SliceChunkReader<'a> {
    label: Option<&'static str>,
    pos: usize,
    line_offset: usize,
    slice: &'a [u8],
    chunk_size: usize,
}

impl<'a> Iterator for SliceChunkReader<'a> {
    type Item = FastqReader<'a, SliceFastqChunkSource<'a>>;
    fn next(&mut self) -> Option<Self::Item> {
        self.chunk_reader()
    }
}

fn slice_chunk<'a>(
    slice: &'a [u8],
    pos: usize,
    chunk_size: usize,
) -> (&'a [u8], Vec<usize>) {
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

    pub fn chunk_reader(
        &mut self,
    ) -> Option<FastqReader<'a, SliceFastqChunkSource<'a>>> {
        // If we've reached or passed the end of the input buffer, stop reading
        if self.pos >= self.slice.len() {
            return None;
        }
        let (mut chunk, mut newlines) =
            slice_chunk(self.slice, self.pos, self.chunk_size);

        // If the chunk contains fewer than 4 lines, the FASTQ record is incomplete
        // so we expand the buffer and try again with a larger chunk
        if newlines.len() < 4 && (self.pos + self.chunk_size) < self.slice.len()
        {
            // we always double the buffer size, if there no enough buffer to hold a whole fastq record
            self.chunk_size *= 2;
            return self.chunk_reader();
        }

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        // SAFETY: above code will ensure newlines won't be empty
        let nkeep = newlines.len() - (newlines.len() % 4);
        newlines.truncate(nkeep);
        let endpoint = unsafe { *newlines.last().unwrap_unchecked() };
        if endpoint == chunk.len() {
            chunk = &chunk[.. endpoint]; // we have no endpoint '\n'
        } else {
            chunk = &chunk[..= endpoint];
        }

        self.pos += chunk.len();
        let line_offset = self.line_offset;
        self.line_offset += newlines.len(); // Increment the number of chunk lines
        let source = SliceFastqChunkSource::new(chunk, newlines);
        let mut out = FastqReader::with_offset(line_offset, source);
        if let Some(label) = self.label {
            out.set_label(label);
        }
        return Some(out);
    }
}

pub struct SliceChunkPairedReader<'a, 'b> {
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
    type Item = std::result::Result<
        FastqPairedReader<
            'a,
            'b,
            SliceFastqChunkSource<'a>,
            SliceFastqChunkSource<'b>,
        >,
        FastqParseError,
    >;
    fn next(&mut self) -> Option<Self::Item> {
        self.chunk_reader().transpose()
    }
}

impl<'a, 'b> SliceChunkPairedReader<'a, 'b> {
    #[allow(dead_code)]
    pub fn new(slice1: &'a [u8], slice2: &'b [u8]) -> Self {
        Self::with_capacity(8 * 1024, slice1, slice2)
    }

    pub fn with_capacity(
        capacity: usize,
        slice1: &'a [u8],
        slice2: &'b [u8],
    ) -> Self {
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

    pub fn chunk_reader(
        &mut self,
    ) -> Result<
        Option<
            FastqPairedReader<
                'a,
                'b,
                SliceFastqChunkSource<'a>,
                SliceFastqChunkSource<'b>,
            >,
        >,
        FastqParseError,
    > {
        let end1 = self.pos1 + self.chunk_size;
        let end2 = self.pos2 + self.chunk_size;

        // Handle EOF logic for paired files
        let eof1 = self.pos1 >= self.slice1.len();
        let eof2 = self.pos2 >= self.slice2.len();

        match (eof1, eof2) {
            (true, true) => return Ok(None),
            (true, false) => {
                if self.slice2[end2 ..].iter().all(|b| b.is_ascii_whitespace())
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
                if self.slice1[end1 ..].iter().all(|b| b.is_ascii_whitespace())
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
            (false, false) => {}
        }

        let (mut chunk1, mut newlines1) =
            slice_chunk(self.slice1, self.pos1, self.chunk_size);

        let (mut chunk2, mut newlines2) =
            slice_chunk(self.slice2, self.pos1, self.chunk_size);

        // Truncate to a multiple of 4 lines to ensure only complete records are included
        let count = usize::min(newlines1.len(), newlines2.len());
        let nkeep = count - (count % 4);
        if nkeep < 4 && end1 < self.slice1.len() && end2 < self.slice2.len() {
            self.chunk_size += self.chunk_size;
            return self.chunk_reader();
        }

        newlines1.truncate(nkeep);
        newlines2.truncate(nkeep);
        let endpoint1 = unsafe { *newlines1.last().unwrap_unchecked() };
        if endpoint1 == chunk1.len() {
            chunk1 = &chunk1[.. endpoint1]; // we have no endpoint '\n'
        } else {
            chunk1 = &chunk1[..= endpoint1];
        }
        let endpoint2 = unsafe { *newlines2.last().unwrap_unchecked() };
        if endpoint2 == chunk2.len() {
            chunk2 = &chunk2[.. endpoint2]; // we have no endpoint '\n'
        } else {
            chunk2 = &chunk2[..= endpoint2];
        }

        self.pos1 += chunk1.len();
        self.pos2 += chunk2.len();

        let line_offset1 = self.line_offset1;
        let line_offset2 = self.line_offset2;

        self.line_offset1 += newlines1.len();
        self.line_offset2 += newlines2.len();

        let mut reader1 = FastqReader::with_offset(
            line_offset1,
            SliceFastqChunkSource::new(chunk1, newlines1),
        );
        let mut reader2 = FastqReader::with_offset(
            line_offset2,
            SliceFastqChunkSource::new(chunk2, newlines2),
        );

        if let Some(label1) = self.label1 {
            reader1.set_label(label1);
        }
        if let Some(label2) = self.label2 {
            reader2.set_label(label2);
        }
        Ok(Some(FastqPairedReader::new(reader1, reader2)))
    }
}

#[derive(Debug)]
pub struct SliceFastqChunkSource<'a> {
    pos: usize,
    newlines: Vec<usize>,
    newlines_pos: usize,
    chunk: &'a [u8],
}

impl<'a> SliceFastqChunkSource<'a> {
    fn new(chunk: &'a [u8], newlines: Vec<usize>) -> Self {
        Self {
            pos: 0,
            newlines,
            newlines_pos: 0,
            chunk,
        }
    }
}

impl<'a> FastqSource<'a> for SliceFastqChunkSource<'a> {
    type Record = FastaRecord<&'a [u8]>;
    fn read_line(&mut self) -> Option<&'a [u8]> {
        if self.newlines_pos >= self.newlines.len() {
            return None;
        }
        let end = unsafe { *self.newlines.get_unchecked(self.newlines_pos) };
        let out = unsafe { self.chunk.get_unchecked(self.pos .. end) };
        self.newlines_pos += 1;
        self.pos = end + 1;
        Some(out)
    }
    fn build_record(
        &mut self,
        id: &'a [u8],
        desc: Option<&'a [u8]>,
        seq: &'a [u8],
        _: &'a [u8],
        _: &'a [u8],
    ) -> Self::Record {
        FastaRecord::new(id, desc, seq)
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
                assert_eq!(r1.id[5..], r2.id[5..]); // match SEQ_IDn
                total_records += 1;
            }
        }
        assert_eq!(total_records, 4);
    }
}
