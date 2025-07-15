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
        label: Option<&'static str>,
        record: String,
        pos: usize,
    },
    UnequalLength {
        label: Option<&'static str>,
        record: String,
        seq: usize,
        qual: usize,
        pos: usize,
    },
    InvalidSep {
        label: Option<&'static str>,
        record: String,
        pos: usize,
    },
    IncompleteRecord {
        label: Option<&'static str>,
        record: String,
        pos: usize,
    },
    FastqPairError {
        read1_label: Option<&'static str>,
        read2_label: Option<&'static str>,
        read1_id: String,
        read2_id: String,
        read1_pos: usize,
        read2_pos: usize,
    },
    OutOfSync {
        eof_label: Option<&'static str>,
        continueal_label: Option<&'static str>,
        eof_pos: usize,
        continueal_pos: usize,
    },
}

impl fmt::Display for FastqParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn label_line(label: Option<&str>, pos: usize) -> String {
            match label {
                Some(label) => format!("in {} (line: {})", label, pos),
                None => format!("(line: {})", pos),
            }
        }

        match self {
            FastqParseError::IncompleteRecord { label, record, pos } => {
                write!(
                    f,
                    "FASTQ parse error {}: incomplete record\n{}",
                    label_line(*label, *pos),
                    record
                )
            }
            FastqParseError::InvalidHead { label, record, pos } => {
                write!(
                    f,
                    "FASTQ parse error {}: expected '@' at record start\n{}",
                    label_line(*label, *pos),
                    record
                )
            }
            FastqParseError::UnequalLength {
                label,
                record,
                seq,
                qual,
                pos,
            } => {
                write!(
                    f,
                    "FASTQ parse error {}: sequence and quality lengths do not match ({} vs {})\n{}",
                    label_line(*label, *pos),
                    seq, qual,
                    record
                )
            }
            FastqParseError::InvalidSep { label, record, pos } => {
                write!(
                    f,
                    "FASTQ parse error {}: expected '+' at separator line start\n{}",
                    label_line(*label, *pos),
                    record
                )
            }
            FastqParseError::FastqPairError {
                read1_label,
                read2_label,
                read1_id,
                read2_id,
                read1_pos,
                read2_pos,
            } => {
                let r1 = read1_label.unwrap_or("read1");
                let r2 = read2_label.unwrap_or("read2");
                write!(
                    f,
                    "FASTQ pairing error: sequence IDs do not match\n  {} (line: {}): {}\n  {} (line: {}): {}",
                    r1, read1_pos, read1_id,
                    r2, read2_pos, read2_id
                )
            }
            FastqParseError::OutOfSync {
                eof_label,
                continueal_label,
                eof_pos,
                continueal_pos,
            } => {
                let r1 = eof_label.unwrap_or("read1");
                let r2 = continueal_label.unwrap_or("read2");
                write!(
                    f,
                    "FASTQ pairing error: reached end of {} (line: {}) while {} (line: {}) still contains additional lines",
                    r1, eof_pos, r2, continueal_pos
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
    fn test_invalid_head_with_label() {
        let err = FastqParseError::InvalidHead {
            label: Some("read1"),
            record: "head: @SEQ_ID".into(),
            pos: 42,
        };
        let msg = format!("{}", err);
        assert!(
            msg.contains("in read1 (line: 42)"),
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
            label: None,
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
    fn test_fastq_pair_error_with_labels() {
        let err = FastqParseError::FastqPairError {
            read1_label: Some("R1"),
            read2_label: Some("R2"),
            read1_id: "SEQ1".into(),
            read2_id: "SEQ2".into(),
            read1_pos: 3,
            read2_pos: 3,
        };
        let msg = format!("{}", err);
        assert!(
            msg.contains("sequence IDs do not match"),
            "Should mention pairing mismatch"
        );
        assert!(msg.contains("R1 (line: 3): SEQ1"));
        assert!(msg.contains("R2 (line: 3): SEQ2"));
    }

    #[test]
    fn test_out_of_sync_default_labels() {
        let err = FastqParseError::OutOfSync {
            eof_label: None,
            continueal_label: None,
            eof_pos: 50,
            continueal_pos: 51,
        };
        let msg = format!("{}", err);
        assert!(msg.contains("read1 (line: 50)"), "Should default to read1");
        assert!(msg.contains("read2 (line: 51)"), "Should default to read2");
    }
}
