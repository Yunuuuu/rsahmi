use std::io::{BufRead, Read, Write};

use indicatif::ProgressBar;

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

impl<R: BufRead> BufRead for ProgressBarReader<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.reader.fill_buf()
    }
    fn consume(&mut self, amount: usize) {
        self.reader.consume(amount);
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
        let n = self.writer.write(buf)?;
        self.bar.inc(n as u64);
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
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
