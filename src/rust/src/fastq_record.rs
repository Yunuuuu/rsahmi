use std::io::Write;

#[derive(Debug)]
pub(crate) struct FastqRecord<T> {
    pub(crate) id: T,
    pub(crate) desc: Option<T>,
    pub(crate) seq: T,
    pub(crate) sep: T,
    pub(crate) qual: T,
}

impl<T> FastqRecord<T> {
    #[allow(dead_code)]
    pub(crate) fn new(id: T, desc: Option<T>, seq: T, sep: T, qual: T) -> Self {
        Self {
            id,
            desc,
            seq,
            sep,
            qual,
        }
    }
}

impl<T: AsRef<[u8]>> FastqRecord<T> {
    #[allow(dead_code)]
    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.as_vec())
    }

    pub(crate) fn bytes_size(&self) -> usize {
        self.id.as_ref().len()
            // extra one for space between id and description
            + self.desc.as_ref().map(|d| d.as_ref().len() + 1).unwrap_or(0) // ' '
            + self.seq.as_ref().len()
            + self.sep.as_ref().len()
            + self.qual.as_ref().len()
            + 5 // '@' and 4 * '\n'
    }

    /// Efficiently appends the FASTQ record to the provided Vec<u8>
    pub(crate) fn extend(&self, buf: &mut Vec<u8>) {
        buf.push(b'@');
        buf.extend_from_slice(self.id.as_ref());

        if let Some(desc) = &self.desc {
            buf.push(b' ');
            buf.extend_from_slice(desc.as_ref());
        }

        buf.push(b'\n');
        buf.extend_from_slice(self.seq.as_ref());
        buf.push(b'\n');
        buf.extend_from_slice(self.sep.as_ref());
        buf.push(b'\n');
        buf.extend_from_slice(self.qual.as_ref());
        buf.push(b'\n');
    }

    pub(crate) fn as_vec(&self) -> Vec<u8> {
        let id = self.id.as_ref();
        let desc = self.desc.as_ref().map(|d| d.as_ref());
        let seq = self.seq.as_ref();
        let sep = self.sep.as_ref();
        let qual = self.qual.as_ref();
        let mut buffer = Vec::with_capacity(self.bytes_size());
        buffer.push(b'@');
        buffer.extend_from_slice(id);

        if let Some(desc) = &desc {
            buffer.push(b' ');
            buffer.extend_from_slice(desc);
        }

        buffer.push(b'\n');
        buffer.extend_from_slice(seq);
        buffer.push(b'\n');
        buffer.extend_from_slice(sep);
        buffer.push(b'\n');
        buffer.extend_from_slice(qual);
        buffer.push(b'\n');
        buffer
    }

    #[allow(dead_code)]
    pub(crate) fn write_buf(&self, buf: &mut [u8]) -> std::io::Result<usize> {
        let mut pos = 0;

        // Write '@' and ID
        buf[pos] = b'@';
        pos += 1;

        let id = self.id.as_ref();
        if pos + id.len() > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + id.len()].copy_from_slice(id);
        pos += id.len();

        // Optional description
        if let Some(desc) = &self.desc {
            if pos + 1 > buf.len() {
                return Err(std::io::ErrorKind::WriteZero.into());
            }
            buf[pos] = b' ';
            pos += 1;

            let desc = desc.as_ref();
            if pos + desc.len() > buf.len() {
                return Err(std::io::ErrorKind::WriteZero.into());
            }
            buf[pos .. pos + desc.len()].copy_from_slice(desc);
            pos += desc.len();
        }

        // Line break
        if pos + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos] = b'\n';
        pos += 1;

        // Sequence
        let seq = self.seq.as_ref();
        if pos + seq.len() + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + seq.len()].copy_from_slice(seq);
        pos += seq.len();
        buf[pos] = b'\n';
        pos += 1;

        // Separator
        let sep = self.sep.as_ref();
        if pos + sep.len() + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + sep.len()].copy_from_slice(sep);
        pos += sep.len();
        buf[pos] = b'\n';
        pos += 1;

        // Quality
        let qual = self.qual.as_ref();
        if pos + qual.len() + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + qual.len()].copy_from_slice(qual);
        pos += qual.len();
        buf[pos] = b'\n';
        pos += 1;

        Ok(pos)
    }

    #[allow(dead_code)]
    pub(crate) fn as_ref(&self) -> FastqRecord<&[u8]> {
        FastqRecord {
            id: self.id.as_ref(),
            desc: self.desc.as_ref().map(|d| d.as_ref()),
            seq: self.seq.as_ref(),
            sep: self.sep.as_ref(),
            qual: self.qual.as_ref(),
        }
    }
}
use std::error::Error;
use std::fmt;

/// FASTQ parsing error
#[derive(Debug)]
pub enum FastqParseError {
    InvalidHead {
        record: String,
        pos: usize,
    },
    UnequalLength {
        record: String,
        seq: usize,
        qual: usize,
        pos: usize,
    },
    InvalidSep {
        record: String,
        pos: usize,
    },
    IncompleteRecord {
        record: String,
        pos: usize,
    },
    FastqPairError {
        read1_id: String,
        read2_id: String,
        read1_pos: Option<usize>,
        read2_pos: Option<usize>,
    },
}

impl fmt::Display for FastqParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FastqParseError::IncompleteRecord { record, pos } => {
                write!(
                    f,
                    "FASTQ parse error (line: {}): incomplete record\n{}",
                    pos, record
                )
            }
            FastqParseError::InvalidHead { record, pos } => {
                write!(
                    f,
                    "FASTQ parse error (line: {}): expected '@' at record start\n{}",
                    pos, record
                )
            }
            FastqParseError::UnequalLength {
                record,
                seq,
                qual,
                pos,
            } => {
                write!(
                    f,
                    "FASTQ parse error (line: {}): sequence and quality lengths do not match ({} vs {})\n{}",
                    pos,
                    seq, qual,
                    record
                )
            }
            FastqParseError::InvalidSep { record, pos } => {
                write!(
                    f,
                    "FASTQ parse error (line: {}): expected '+' at separator line start\n{}",
                    pos, record
                )
            }
            FastqParseError::FastqPairError {
                read1_id,
                read2_id,
                read1_pos,
                read2_pos,
            } => {
                write!(
                    f,
                    "FASTQ pairing error: sequence IDs do not match\n  record1 ID{}: {}\n  record2 ID{}: {}",
                    match read1_pos {
                        Some(pos) => format!(" (line: {})", pos),
                        None => "".to_string(),
                    },
                    read1_id,
                    match read2_pos {
                        Some(pos) => format!(" (line: {})", pos),
                        None => "".to_string(),
                    },
                    read2_id
                )
            }
        }
    }
}

impl Error for FastqParseError {}

#[cfg(test)]
mod test_record {
    use std::io::Cursor;

    use super::*; // Import FastqRecord

    #[test]
    fn test_fastq_record_write_with_description() {
        let record = FastqRecord::new(
            b"SEQ_ID".as_ref(),
            Some(b"desc".as_ref()),
            b"ACGTACGT".as_ref(),
            b"+".as_ref(),
            b"IIIIIIII".as_ref(),
        );

        let mut output = Cursor::new(Vec::new());
        record.write(&mut output).expect("Write failed");

        let expected = b"@SEQ_ID desc\nACGTACGT\n+\nIIIIIIII\n";
        assert_eq!(output.into_inner(), expected);
    }

    #[test]
    fn test_fastq_record_write_without_description() {
        let record = FastqRecord::new(
            b"SEQ_ID".as_ref(),
            None,
            b"ACGTACGT".as_ref(),
            b"+".as_ref(),
            b"IIIIIIII".as_ref(),
        );

        let mut output = Cursor::new(Vec::new());
        record.write(&mut output).expect("Write failed");

        let expected = b"@SEQ_ID\nACGTACGT\n+\nIIIIIIII\n";
        assert_eq!(output.into_inner(), expected);
    }
}

#[cfg(test)]
mod test_error {
    use super::*;

    #[test]
    fn test_invalid_head() {
        let err = FastqParseError::InvalidHead {
            record: "head: @SEQ_ID".into(),
            pos: 42,
        };
        let msg = format!("{}", err);
        assert!(
            msg.contains("(line: 42)"),
            "Message should contain label and line"
        );
        assert!(
            msg.contains("expected '@' at record start"),
            "Message should mention invalid head"
        );
    }

    #[test]
    fn test_unequal_length_without_label() {
        let err = FastqParseError::UnequalLength {
            record: "SEQ\nQUAL".into(),
            seq: 10,
            qual: 8,
            pos: 100,
        };
        let msg = format!("{}", err);
        assert!(
            msg.contains("(line: 100)"),
            "Should fall back to just line number"
        );
        assert!(
            msg.contains("sequence and quality lengths do not match (10 vs 8)"),
            "Should indicate length mismatch"
        );
    }

    #[test]
    fn test_fastq_pair_error() {
        let err = FastqParseError::FastqPairError {
            read1_id: "SEQ1".into(),
            read2_id: "SEQ2".into(),
            read1_pos: Some(3),
            read2_pos: Some(3),
        };
        let msg = format!("{}", err);
        assert!(
            msg.contains("sequence IDs do not match"),
            "Should mention pairing mismatch"
        );
        assert!(msg.contains("(line: 3): SEQ1"));
        assert!(msg.contains("(line: 3): SEQ2"));
    }
}
