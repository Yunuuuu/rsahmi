use std::io::Read;

use anyhow::Result;
use bytes::{Bytes, BytesMut};
use memchr::memchr2;

use crate::fastq_record::FastqParseError;
use crate::fastq_record::FastqRecord;
use crate::reader0::*;

pub(crate) struct FastqReader<R> {
    reader: LineReader<R>,
}

impl<R: Read> FastqReader<R> {
    #[allow(dead_code)]
    pub(crate) fn new(reader: R) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub(crate) fn with_capacity(capacity: usize, reader: R) -> Self {
        Self {
            reader: LineReader::with_capacity(capacity, reader),
        }
    }

    pub(crate) fn offset(&self) -> usize {
        self.reader.offset()
    }

    #[inline]
    fn read_line(&mut self) -> std::io::Result<Option<BytesMut>> {
        self.reader.read_line()
    }

    #[inline]
    pub(crate) fn read_record(&mut self) -> Result<Option<FastqRecord<Bytes>>> {
        let mut header;
        loop {
            if let Some(line) = self.read_line()? {
                if line.iter().all(|b| b.is_ascii_whitespace()) {
                    continue;
                } else {
                    header = line;
                    break;
                }
            } else {
                return Ok(None);
            }
        }

        // SAFETY: we must ensure line is not empty, this is ensured by the caller function
        if header.is_empty() || unsafe { *header.get_unchecked(0) } != b'@' {
            Err(FastqParseError::InvalidHead {
                label: None,
                record: format!("{}", String::from_utf8_lossy(&header)),
                pos: self.offset(),
            })?;
        }
        let _ = header.split_to(1); // remove the '@' from the start of the sequence ID
        let id;
        let desc;
        if let Some(line_pos) = memchr2(b' ', b'\t', &header) {
            id = header.split_to(line_pos).freeze();
            let _ = header.split_to(1); // remove the blankspace
                                        // check if description exits
            if header.is_empty() {
                desc = None
            } else {
                desc = Some(header.freeze());
            }
        } else {
            id = header.freeze();
            desc = None;
        }

        // 2nd line (sequence) must exist. Otherwise, incomplete record.
        let seq = if let Some(line) = self.read_line()? {
            Ok(line.freeze())
        } else {
            Err(FastqParseError::IncompleteRecord {
                label: None,
                record: format!(
                    "{}{}",
                    String::from_utf8_lossy(&id),
                    String::from_utf8_lossy(
                        desc.as_ref()
                            .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                    )
                ),
                pos: self.offset(),
            })
        }?;

        // 3rd line (separator). Must exist.
        let sep = if let Some(line) = self.read_line()? {
            // Separator: begins with a '+' character and is optionally followed by the same sequence identifier
            if line.is_empty() || unsafe { *line.get_unchecked(0) } != b'+' {
                Err(FastqParseError::InvalidSep {
                    label: None,
                    record: format!(
                        "{}{}\n{}\n{}",
                        String::from_utf8_lossy(&id),
                        String::from_utf8_lossy(
                            desc.as_ref()
                                .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                        ),
                        String::from_utf8_lossy(&seq),
                        String::from_utf8_lossy(&line)
                    ),
                    pos: self.offset(),
                })
            } else {
                Ok(line.freeze())
            }
        } else {
            Err(FastqParseError::IncompleteRecord {
                label: None,
                record: format!(
                    "{}{}\n{}",
                    String::from_utf8_lossy(&id),
                    String::from_utf8_lossy(
                        desc.as_ref()
                            .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                    ),
                    String::from_utf8_lossy(&seq)
                ),
                pos: self.offset(),
            })
        }?;

        // 4th line (quality). Must exist.
        let qual = if let Some(line) = self.read_line()? {
            if seq.len() != line.len() {
                Err(FastqParseError::UnequalLength {
                    label: None,
                    seq: seq.len(),
                    qual: line.len(),
                    record: format!(
                        "{}{}\n{}\n{}\n{}",
                        String::from_utf8_lossy(&id),
                        String::from_utf8_lossy(
                            desc.as_ref()
                                .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                        ),
                        String::from_utf8_lossy(&seq),
                        String::from_utf8_lossy(&sep),
                        String::from_utf8_lossy(&line)
                    ),
                    pos: self.offset(),
                })
            } else {
                Ok(line.freeze())
            }
        } else {
            Err(FastqParseError::IncompleteRecord {
                label: None,
                record: format!(
                    "{}\n{}\n{}",
                    String::from_utf8_lossy(
                        desc.as_ref()
                            .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                    ),
                    String::from_utf8_lossy(&seq),
                    String::from_utf8_lossy(&sep),
                ),
                pos: self.offset(),
            })
        }?;
        Ok(Some(FastqRecord::new(id, desc, seq, sep, qual)))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn create_reader(data: &str) -> FastqReader<Cursor<&[u8]>> {
        let reader = Cursor::new(data.as_bytes());
        FastqReader::new(reader)
    }

    #[test]
    fn test_read_valid_record() -> Result<()> {
        let fastq_data = "@seq1 description\nATGC\n+\n!!!!\n@seq2\nGCGT\n+\n$$$$\n";

        let mut reader = create_reader(fastq_data);

        let record = reader.read_record()?.expect("Should have a record");

        // Check that the FASTQ fields match the expected values
        assert_eq!(record.id.as_ref(), b"seq1");
        assert_eq!(record.desc.unwrap().as_ref(), b"description");
        assert_eq!(record.seq.as_ref(), b"ATGC");
        assert_eq!(record.sep.as_ref(), b"+");
        assert_eq!(record.qual.as_ref(), b"!!!!");

        Ok(())
    }

    #[test]
    fn test_invalid_header() -> Result<()> {
        let fastq_data = "seq1 description\nATGC\n+\n!!!!\n";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_record();

        assert!(result.is_err()); // Expect error: Invalid header
        Ok(())
    }

    #[test]
    fn test_incomplete_record() -> Result<()> {
        let fastq_data = "@seq1 description\nATGC\n+\n";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_record();

        assert!(result.is_err()); // Expect error: Incomplete record (missing quality)
        Ok(())
    }

    #[test]
    fn test_unmatched_seq_qual_length() -> Result<()> {
        let fastq_data = "@seq1 description\nATGC\n+\n!!\n";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_record();

        assert!(result.is_err()); // Expect error: Unequal sequence and quality lengths
        Ok(())
    }

    #[test]
    fn test_edge_case_empty_data() -> Result<()> {
        let fastq_data = "";

        let mut reader = create_reader(fastq_data);

        let result = reader.read_record();

        assert!(result.is_ok()); // Expect error: EOF
        assert!(result.unwrap().is_none());
        Ok(())
    }
}
