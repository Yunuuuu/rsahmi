use std::fs::File;
use std::io::BufReader;
use std::io::{BufRead, Read, Write};
use std::path::Path;

use anyhow::{anyhow, Result};
use bytes::{Bytes, BytesMut};
use indicatif::ProgressBar;
use isal::read::GzipDecoder;
use libdeflater::Compressor;
use memchr::{memchr, memchr2};

use crate::parser::fastq::FastqParseError;
use crate::parser::fastq::FastqRecord;
use crate::reader0::*;

pub(crate) struct FastqReader<R> {
    reader: R,
    offset: usize,
    leftover: Option<BytesMut>,
    label: Option<&'static str>,
}

pub(crate) fn gz_compressed(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .map_or(false, |s| s.eq_ignore_ascii_case("gz"))
}

pub(crate) fn fastq_writer<P: AsRef<Path> + ?Sized>(
    file: &P,
    progress_bar: Option<ProgressBar>,
) -> Result<Box<dyn Write>> {
    let path: &Path = file.as_ref();
    let file = File::create(path)
        .map_err(|e| anyhow!("Failed to create output file {}: {}", path.display(), e))?;
    let writer: Box<dyn Write>;
    if let Some(bar) = progress_bar {
        writer = Box::new(ProgressBarWriter::new(file, bar));
    } else {
        writer = Box::new(file);
    }
    Ok(writer)
}

pub(crate) fn fastq_pack(bytes: &[u8], compressor: &mut Compressor, gzip: bool) -> Result<Vec<u8>> {
    let mut pack: Vec<u8>;
    if gzip {
        let pack_size = compressor.gzip_compress_bound(bytes.len());
        pack = Vec::with_capacity(pack_size);
        pack.resize(pack_size, 0);
        let size = compressor.gzip_compress(bytes, &mut pack)?;
        pack.truncate(size);
    } else {
        pack = Vec::with_capacity(bytes.len());
        pack.extend_from_slice(bytes);
    }
    Ok(pack)
}

pub(crate) fn fastq_reader<P: AsRef<Path> + ?Sized>(
    file: &P,
    buffer_size: usize,
    progress_bar: Option<ProgressBar>,
) -> Result<Box<dyn Read>> {
    let path: &Path = file.as_ref();
    let file =
        File::open(path).map_err(|e| anyhow!("Failed to open file {}: {}", path.display(), e))?;
    let reader: Box<dyn Read>;
    if gz_compressed(path) {
        if let Some(bar) = progress_bar {
            reader = Box::new(GzipDecoder::new(BufReader::with_capacity(
                buffer_size,
                ProgressBarReader::new(file, bar),
            )));
        } else {
            reader = Box::new(GzipDecoder::new(BufReader::with_capacity(
                buffer_size,
                file,
            )));
        }
    } else {
        if let Some(bar) = progress_bar {
            reader = Box::new(ProgressBarReader::new(file, bar));
        } else {
            reader = Box::new(file);
        }
    }
    Ok(reader)
}

impl<R: Read> FastqReader<BufReader<R>> {
    pub(crate) fn new(reader: BufReader<R>) -> Self {
        Self {
            reader,
            offset: 0,
            leftover: None,
            label: None,
        }
    }

    #[inline]
    fn read_line(&mut self) -> Result<Option<BytesMut>> {
        loop {
            let buffer = self.reader.fill_buf()?;
            if buffer.is_empty() {
                return Ok(self.leftover.take().filter(|b| !b.is_empty()));
            }

            if let Some(pos) = memchr(b'\n', buffer) {
                // Fast path: newline found
                let end = if pos > 0 && buffer[pos - 1] == b'\r' {
                    pos - 1
                } else {
                    pos
                };

                let line = if let Some(mut leftover) = self.leftover.take() {
                    leftover.extend_from_slice(&buffer[.. end]);
                    leftover
                } else {
                    // Directly build from slice without heap copying if possible
                    BytesMut::from(&buffer[.. end])
                };

                self.reader.consume(pos + 1);
                return Ok(Some(line));
            }

            // No newline: accumulate buffer and continue
            let len = buffer.len();
            if self.leftover.is_none() {
                self.leftover = Some(BytesMut::with_capacity(len.max(64)));
            }
            self.leftover.as_mut().unwrap().extend_from_slice(buffer);
            self.reader.consume(len);
        }
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
                label: self.label,
                record: format!("{}", String::from_utf8_lossy(&header)),
                pos: self.offset,
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
                label: self.label,
                record: format!(
                    "{}{}",
                    String::from_utf8_lossy(&id),
                    String::from_utf8_lossy(
                        desc.as_ref()
                            .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                    )
                ),
                pos: self.offset,
            })
        }?;

        // 3rd line (separator). Must exist.
        let sep = if let Some(line) = self.read_line()? {
            // Separator: begins with a '+' character and is optionally followed by the same sequence identifier
            if line.is_empty() || unsafe { *line.get_unchecked(0) } != b'+' {
                Err(FastqParseError::InvalidSep {
                    label: self.label,
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
                    pos: self.offset,
                })
            } else {
                Ok(line.freeze())
            }
        } else {
            Err(FastqParseError::IncompleteRecord {
                label: self.label,
                record: format!(
                    "{}{}\n{}",
                    String::from_utf8_lossy(&id),
                    String::from_utf8_lossy(
                        desc.as_ref()
                            .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                    ),
                    String::from_utf8_lossy(&seq)
                ),
                pos: self.offset,
            })
        }?;

        // 4th line (quality). Must exist.
        let qual = if let Some(line) = self.read_line()? {
            if seq.len() != line.len() {
                Err(FastqParseError::UnequalLength {
                    label: self.label,
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
                    pos: self.offset,
                })
            } else {
                Ok(line.freeze())
            }
        } else {
            Err(FastqParseError::IncompleteRecord {
                label: self.label,
                record: format!(
                    "{}\n{}\n{}",
                    String::from_utf8_lossy(
                        desc.as_ref()
                            .map_or_else(|| -> &[u8] { b"" }, |d| -> &[u8] { &d })
                    ),
                    String::from_utf8_lossy(&seq),
                    String::from_utf8_lossy(&sep),
                ),
                pos: self.offset,
            })
        }?;
        Ok(Some(FastqRecord::new(id, desc, seq, sep, qual)))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn create_reader(data: &str) -> FastqReader<BufReader<Cursor<&[u8]>>> {
        let reader = Cursor::new(data.as_bytes());
        let bytes_reader = BufReader::new(reader);
        FastqReader {
            reader: bytes_reader,
            offset: 0,
            leftover: None,
            label: None,
        }
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
