use memchr::memchr2;

use super::FastqContainer;
use super::FastqParseError;

pub trait FastqSource<'a>: Sized {
    type Record;

    fn read_line(&mut self) -> Option<&'a [u8]>;

    fn build_record(
        &mut self,
        id: &'a [u8],
        desc: Option<&'a [u8]>,
        seq: &'a [u8],
        sep: &'a [u8],
        qual: &'a [u8],
    ) -> Self::Record;

    #[allow(dead_code)]
    fn into_reader(self) -> FastqReader<'a, Self> {
        FastqReader::new(self)
    }
}

#[derive(Debug)]
pub struct FastqReader<'a, S>
where
    S: FastqSource<'a>,
{
    label: Option<&'static str>,
    source: S,
    offset: usize,
    container: FastqContainer<'a>,
}

impl<'a, S> FastqReader<'a, S>
where
    S: FastqSource<'a>,
{
    pub fn new(source: S) -> Self {
        Self::with_offset(0, source)
    }

    pub fn with_offset(offset: usize, source: S) -> Self {
        Self {
            label: None,
            source,
            offset,
            container: FastqContainer::default(),
        }
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

    #[allow(dead_code)]
    pub fn into_reader(self) -> S {
        self.source
    }

    pub fn read_record(
        &mut self,
    ) -> Result<Option<S::Record>, FastqParseError> {
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

    pub fn id(&self) -> Option<&'a [u8]> {
        self.container.id()
    }

    pub fn desc(&self) -> Option<Option<&'a [u8]>> {
        self.container.desc()
    }

    pub fn seq(&self) -> Option<&'a [u8]> {
        self.container.seq()
    }

    pub fn sep(&self) -> Option<&'a [u8]> {
        self.container.sep()
    }

    pub fn qual(&self) -> Option<&'a [u8]> {
        self.container.qual()
    }

    // --- internal method, must call for cautious, they use unsafe code;
    fn parse_head(&mut self, line: &'a [u8]) -> Result<(), FastqParseError> {
        // SAFETY: we must ensure line is not empty, this is ensured by the caller function
        if unsafe { *line.get_unchecked(0) } != b'@' {
            return Err(FastqParseError::InvalidHead {
                label: self.label,
                record: format!("head: {}", String::from_utf8_lossy(line)),
                pos: self.offset,
            });
        }
        let id;
        let desc;
        if let Some(line_pos) = memchr2(b' ', b'\t', line) {
            // remove the '@' from the start of the sequence ID
            id = &line[1 .. line_pos];
            // check if description exits
            if line_pos + 1 == line.len() {
                desc = None
            } else {
                desc = Some(&line[line_pos + 1 ..]);
            }
        } else {
            // remove the '@' from the start of the sequence ID
            id = &line[1 ..];
            desc = None;
        }
        self.container.set_id(id);
        self.container.set_desc(desc);
        Ok(())
    }

    fn parse_seq(&mut self, line: &'a [u8]) -> Result<(), FastqParseError> {
        self.container.set_seq(line);
        Ok(())
    }

    fn parse_sep(&mut self, line: &'a [u8]) -> Result<(), FastqParseError> {
        // Separator: begins with a '+' character and is optionally followed by the same sequence identifier
        if line.is_empty() || unsafe { *line.get_unchecked(0) } != b'+' {
            return Err(FastqParseError::InvalidSep {
                label: self.label,
                record: format!(
                    "{}\nseparator: {}",
                    self.container.to_string(),
                    String::from_utf8_lossy(line)
                ),
                pos: self.offset,
            });
        }
        self.container.set_sep(line);
        Ok(())
    }

    fn parse_qual(&mut self, line: &'a [u8]) -> Result<(), FastqParseError> {
        let seq = unsafe { self.seq().unwrap_unchecked() };
        if seq.len() != line.len() {
            return Err(FastqParseError::UnequalLength {
                label: self.label,
                seq: seq.len(),
                qual: line.len(),
                record: format!(
                    "{}\nquality: {}",
                    self.container.to_string(),
                    String::from_utf8_lossy(line)
                ),
                pos: self.offset,
            });
        }
        self.container.set_qual(line);
        Ok(())
    }

    fn read_head(&mut self) -> Result<Option<()>, FastqParseError> {
        // Try reading the 1st line (head). If EOF, return None.
        if let Some(line) = self.source.read_line() {
            // skip the empty line
            if line.is_empty() || line.iter().all(|b| b.is_ascii_whitespace()) {
                self.offset += 1;
                return self.read_head();
            }
            self.parse_head(line)?;
            self.offset += 1;
            Ok(Some(()))
        } else {
            return Ok(None);
        }
    }

    fn read_tail(&mut self) -> Result<(), FastqParseError> {
        // 2nd line (sequence) must exist. Otherwise, incomplete record.
        if let Some(line) = self.source.read_line() {
            self.parse_seq(line)?;
            self.offset += 1;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.label,
                record: self.container.to_string(),
                pos: self.offset,
            });
        }

        // 3rd line (separator). Must exist.
        if let Some(line) = self.source.read_line() {
            self.parse_sep(line)?;
            self.offset += 1;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.label,
                record: self.container.to_string(),
                pos: self.offset,
            });
        }

        // 4th line (quality). Must exist.
        if let Some(line) = self.source.read_line() {
            self.parse_qual(line)?;
            self.offset += 1;
        } else {
            return Err(FastqParseError::IncompleteRecord {
                label: self.label,
                record: self.container.to_string(),
                pos: self.offset,
            });
        }
        Ok(())
    }

    fn build(&mut self) -> S::Record {
        let out = unsafe {
            self.source.build_record(
                self.id().unwrap_unchecked(),
                self.desc().unwrap_unchecked(),
                self.seq().unwrap_unchecked(),
                self.sep().unwrap_unchecked(),
                self.qual().unwrap_unchecked(),
            )
        };
        self.container.reset();
        out
    }
}

pub struct FastqPairedReader<'a, 'b, Reader1, Reader2>
where
    Reader1: FastqSource<'a>,
    Reader2: FastqSource<'b>,
{
    reader1: FastqReader<'a, Reader1>,
    reader2: FastqReader<'b, Reader2>,
}

impl<'a, 'b, Reader1, Reader2> FastqPairedReader<'a, 'b, Reader1, Reader2>
where
    Reader1: FastqSource<'a>,
    Reader2: FastqSource<'b>,
{
    pub fn new(
        reader1: FastqReader<'a, Reader1>,
        reader2: FastqReader<'b, Reader2>,
    ) -> Self {
        Self { reader1, reader2 }
    }

    pub fn read_record(
        &mut self,
    ) -> Result<Option<(Reader1::Record, Reader2::Record)>, FastqParseError>
    {
        match (self.reader1.read_head()?, self.reader2.read_head()?) {
            (Some(()), None) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader2.label,
                    continueal_label: self.reader1.label,
                    eof_pos: self.reader2.offset,
                    continueal_pos: self.reader1.offset,
                });
            }
            (None, Some(())) => {
                return Err(FastqParseError::OutOfSync {
                    eof_label: self.reader1.label,
                    continueal_label: self.reader2.label,
                    eof_pos: self.reader1.offset,
                    continueal_pos: self.reader2.offset,
                });
            }
            (Some(()), Some(())) => {
                let id1 = unsafe { self.reader1.id().unwrap_unchecked() };
                let id2 = unsafe { self.reader2.id().unwrap_unchecked() };
                if id1 != id2 {
                    return Err(FastqParseError::FastqPairError {
                        read1_label: self.reader1.label,
                        read2_label: self.reader2.label,
                        read1_id: String::from_utf8_lossy(id1).into_owned(),
                        read2_id: String::from_utf8_lossy(id2).into_owned(),
                        read1_pos: self.reader1.offset,
                        read2_pos: self.reader2.offset,
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
