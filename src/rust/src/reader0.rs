use std::io::Read;

use bytes::{Bytes, BytesMut};
use indicatif::ProgressBar;
use memchr::memchr;

pub(crate) struct ProgressBarReader<R: Read> {
    bar: ProgressBar,
    reader: R,
}

impl<R: Read> ProgressBarReader<R> {
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

// For dynamic and zero-copy buffer
#[derive(Debug)]
pub(crate) struct BytesReader<R> {
    reader: R,
    buffer_size: usize,
    buffer: BytesMut,
}

impl<R: Read> BytesReader<R> {
    pub(crate) fn new(reader: R) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub(crate) fn with_capacity(capacity: usize, reader: R) -> Self {
        Self {
            reader,
            buffer_size: capacity,
            buffer: BytesMut::with_capacity(capacity),
        }
    }

    // pub(crate) fn label(&self) -> Option<&'static str> {
    //     self.label
    // }
    // pub(crate) fn set_label(&mut self, label: &'static str) {
    //     self.label = Some(label);
    // }
    // #[allow(dead_code)]
    // pub(crate) fn unset_label(&mut self) {
    //     self.label = None;
    // }

    pub(crate) fn into_inner(self) -> R {
        self.reader
    }

    pub(crate) fn borrow_reader(&mut self) -> &mut R {
        self.reader.by_ref()
    }

    pub(crate) fn borrow_buf(&self) -> &BytesMut {
        &self.buffer
    }

    pub(crate) fn borrow_buf_mut(&mut self) -> &mut BytesMut {
        &mut self.buffer
    }

    pub(crate) fn fill_full(&mut self) -> std::io::Result<usize> {
        let leftover_len = self.buffer.len();
        let n_read;
        if leftover_len < self.buffer_size {
            self.buffer.reserve(self.buffer_size - leftover_len);
            unsafe { self.buffer.set_len(self.buffer_size) };
            n_read = self.reader.read(&mut self.buffer[leftover_len ..])?;
            unsafe { self.buffer.set_len(leftover_len + n_read) };
        } else {
            n_read = 0;
        }
        Ok(n_read)
    }

    pub(crate) fn take(&mut self, size: usize) -> std::io::Result<Option<BytesMut>> {
        // if no data left, we just check if there are any data in the underlying reader
        if self.buffer.len() == 0 {
            let n_read = self.fill_full()?;
            if n_read == 0 {
                return Ok(None);
            }
        }

        if size == 0 {
            return Ok(Some(BytesMut::new()));
        }

        // Enough data available, split and return
        if self.buffer.len() >= size {
            return Ok(Some(self.buffer.split_to(size)));
        }

        // No enough data; take what's left
        let mut buf =
            std::mem::replace(&mut self.buffer, BytesMut::with_capacity(self.buffer_size));
        let mut need = size - self.buffer.len();
        buf.reserve(need);

        // Refill until we get enough or EOF
        loop {
            let n = self.fill_full()?;
            if n > 0 {
                if n >= need {
                    // add additional data
                    buf.extend_from_slice(&self.buffer.split_to(need));
                    break;
                } else {
                    buf.extend_from_slice(&self.buffer);
                    self.buffer.clear();
                    need = need - n;
                }
            } else {
                break;
            }
        }
        Ok(Some(buf))
    }

    pub(crate) fn read_until(&mut self, byte: u8) -> std::io::Result<Option<BytesMut>> {
        // if no data left, we just check if there are any data in the underlying reader
        if self.buffer.len() == 0 {
            let n_read = self.fill_full()?;
            if n_read == 0 {
                return Ok(None);
            }
        }

        // Enough data available, split and return
        if let Some(pos) = memchr(byte, &self.buffer) {
            return Ok(Some(self.buffer.split_to(pos + 1)));
        }

        // No enough data; take what's left
        let mut buf =
            std::mem::replace(&mut self.buffer, BytesMut::with_capacity(self.buffer_size));

        // Refill until we get the `byte` or EOF
        loop {
            let n = self.fill_full()?;
            if n > 0 {
                if let Some(pos) = memchr(byte, &self.buffer) {
                    buf.extend_from_slice(&self.buffer.split_to(pos + 1));
                } else {
                    buf.extend_from_slice(&self.buffer);
                    self.buffer.clear();
                }
            } else {
                break;
            }
        }
        Ok(Some(buf))
    }

    pub(crate) fn read_line(&mut self) -> std::io::Result<Option<BytesMut>> {
        self.read_until(b'\n')
    }
}

#[derive(Debug)]
pub(crate) struct BytesDeflateReader<D, R> {
    reader: BytesReader<R>,
    buffer: BytesMut,
    deflater: D,
}

pub(crate) trait Deflate {
    fn deflate(&self, bytes: &[u8]) -> Bytes;
}

impl<D: Deflate, R: Read> BytesDeflateReader<D, R> {
    fn new(deflater: D, reader: R) -> Self {
        Self {
            reader: BytesReader::new(reader),
            buffer: BytesMut::new(),
            deflater,
        }
    }
}
#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::io::{self, Read};

    use bytes::{Bytes, BytesMut};
    use indicatif::ProgressBar;

    use super::{BytesReader, ProgressBarReader};

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

    #[test]
    fn test_bytes_reader_take() {
        let data = get_test_data();
        let cursor = Cursor::new(&data); // Borrow the data here

        let mut reader = BytesReader::new(cursor);
        let chunk_size = 10;
        let expected_data = &data[.. chunk_size];

        // Read the first chunk of data
        if let Some(chunk) = reader.take(chunk_size).unwrap() {
            assert_eq!(chunk, expected_data);
        } else {
            panic!("Failed to read expected chunk of data.");
        }

        // Read the remaining data
        if let Some(chunk) = reader.take(chunk_size).unwrap() {
            assert_eq!(&chunk[..], &data[chunk_size .. chunk_size * 2]);
        } else {
            panic!("Failed to read the remaining chunk of data.");
        }
    }

    #[test]
    fn test_bytes_reader_read_until_newline() {
        let data = b"Hello, world!\nThis is another test.";
        let cursor = Cursor::new(data); // Borrow the data here

        let mut reader = BytesReader::new(cursor);

        // Test that it reads up to the newline character
        if let Some(chunk) = reader.read_line().unwrap() {
            let expected = b"Hello, world!\n";
            assert_eq!(&chunk[..], expected);
        } else {
            panic!("Failed to read up to newline.");
        }
    }

    #[test]
    fn test_bytes_reader_fill_full() {
        let data = get_test_data();
        let cursor = Cursor::new(&data); // Borrow the data here
        let mut reader = BytesReader::new(cursor);

        // Fill the buffer and check if it is filled correctly
        reader.fill_full().unwrap();
        assert_eq!(reader.borrow_buf().len(), data.len());
    }

    #[test]
    fn test_empty_input() {
        let data = Vec::<u8>::new(); // Empty input
        let cursor = Cursor::new(&data); // Borrow the data here
        let mut reader = BytesReader::new(cursor);

        // Try to read from an empty buffer
        let result = reader.read_line().unwrap();
        assert!(result.is_none()); // It should return None (EOF)
    }
}
