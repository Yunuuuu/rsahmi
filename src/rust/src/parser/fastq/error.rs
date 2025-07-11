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
mod tests {
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
