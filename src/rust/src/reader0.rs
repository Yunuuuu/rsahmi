use std::io::{Read, Write};

use bytes::BytesMut;
use indicatif::ProgressBar;
use memchr::memchr;

pub(crate) struct ProgressBarReader<R> {
    bar: ProgressBar,
    reader: R,
}

impl<R> ProgressBarReader<R> {
    pub(crate) fn new(reader: R, bar: ProgressBar) -> Self {
        Self { bar, reader }
    }
}

impl<R: Read> Read for ProgressBarReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let nbytes = self.reader.read(buf)?;
        self.bar.inc(nbytes as u64);
        Ok(nbytes)
    }
}

pub(crate) struct ProgressBarWriter<W> {
    bar: ProgressBar,
    writer: W,
}

impl<W> ProgressBarWriter<W> {
    pub(crate) fn new(writer: W, bar: ProgressBar) -> Self {
        Self { bar, writer }
    }
}

impl<W: Write> Write for ProgressBarWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let nbytes = self.writer.write(buf)?;
        self.bar.inc(nbytes as u64);
        Ok(nbytes)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}

/// LineReader: Efficient zero-copy line-based reader using BytesMut.
///
/// This reader avoids unnecessary heap allocations and copying by:
/// - Reusing a fixed-size buffer (`BytesMut`)
/// - Using `split_to()` to transfer ownership without copying
/// - Accumulating "leftover" when a line spans multiple reads
///
/// Supports CRLF or LF endings and returns each line as a `BytesMut`.
pub(crate) struct LineReader<R> {
    reader: R,                  // Underlying reader (e.g., File)
    offset: usize,              // Line count
    buffer_size: usize,         // buffer capacity
    buffer: Option<BytesMut>,   // Current buffer filled from reader
    leftover: Option<BytesMut>, // Accumulates data when line spans multiple buffers
}

impl<R: Read> LineReader<R> {
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn new(reader: R) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    #[inline]
    pub(crate) fn with_capacity(capacity: usize, reader: R) -> Self {
        Self {
            reader,
            offset: 0,
            buffer: None,
            buffer_size: capacity,
            leftover: None,
        }
    }

    #[inline]
    pub(crate) fn offset(&self) -> usize {
        self.offset
    }

    #[inline]
    pub(crate) fn read_line(&mut self) -> std::io::Result<Option<BytesMut>> {
        loop {
            self.fill_buf()?;
            if let Some(buffer) = self.buffer.as_mut() {
                if let Some(pos) = memchr(b'\n', &buffer) {
                    // Fast path: newline found
                    let mut buf = buffer.split_to(pos + 1);
                    let end = if pos > 0 && buf[pos - 1] == b'\r' {
                        pos - 1
                    } else {
                        pos
                    };
                    let line = if let Some(mut leftover) = self.leftover.take() {
                        leftover.extend_from_slice(&buf[.. end]);
                        leftover
                    } else {
                        // Directly build from slice without heap copying if possible
                        buf.split_to(end)
                    };
                    self.offset += 1;
                    return Ok(Some(line));
                }

                // No newline: accumulate leftover and continue
                if let Some(left) = self.leftover.as_mut() {
                    left.extend_from_slice(&buffer);
                    self.buffer = None
                } else {
                    std::mem::swap(&mut self.buffer, &mut self.leftover);
                }
            } else {
                let left = std::mem::take(&mut self.leftover);
                if left.is_some() {
                    self.offset += 1;
                }
                return Ok(left);
            }
        }
    }

    #[inline]
    fn fill_buf(&mut self) -> std::io::Result<()> {
        if self.buffer.is_none() {
            let mut buffer = BytesMut::with_capacity(self.buffer_size);
            unsafe { buffer.set_len(self.buffer_size) };
            let nbytes = self.reader.read(&mut buffer)?;
            unsafe { buffer.set_len(nbytes) };
            if nbytes > 0 {
                self.buffer = Some(buffer)
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::io::Read;

    use indicatif::ProgressBar;

    use super::ProgressBarReader;

    // Mock input for testing
    fn get_test_data() -> Vec<u8> {
        b"Hello, this is a test input. We will use it to check the readers.".to_vec()
    }

    #[test]
    fn test_progress_bar_reader() {
        let data = get_test_data();
        let cursor = Cursor::new(&data); // Borrow the data here

        // Create progress bar for the reader
        let pb = ProgressBar::new(data.len() as u64);
        let mut reader = ProgressBarReader::new(cursor, pb.clone());

        let mut buffer = vec![0; data.len()];
        let bytes_read = reader.read(&mut buffer).unwrap();

        // Check that we have read the correct number of bytes
        assert_eq!(bytes_read, data.len());

        // Check that the data in the buffer matches the expected input
        assert_eq!(buffer, data);

        // The progress bar should have updated correctly
        assert_eq!(pb.position(), data.len() as u64);
    }
}
