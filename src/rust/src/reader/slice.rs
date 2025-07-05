use indicatif::ProgressBar;
use memchr::{memchr_iter, Memchr};

#[derive(Debug)]
pub(crate) struct SliceReader<'a> {
    slice: &'a [u8],
    pos: usize,
    label: Option<&'static str>,
}

#[derive(Debug)]
pub(crate) struct SliceProgressBarReader<'a> {
    bar: Option<ProgressBar>,
    reader: SliceReader<'a>,
}

impl<'a> SliceProgressBarReader<'a> {
    pub(crate) fn new(slice: &'a [u8]) -> Self {
        Self {
            bar: None,
            reader: SliceReader::new(slice),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn attach_bar(&mut self, bar: ProgressBar) {
        self.bar = Some(bar);
    }

    #[allow(dead_code)]
    pub(crate) fn detach_bar(&mut self) {
        self.bar = None
    }

    pub(crate) fn label(&self) -> Option<&'static str> {
        self.reader.label()
    }

    pub(crate) fn set_label(&mut self, label: &'static str) {
        self.reader.set_label(label)
    }

    #[allow(dead_code)]
    pub(crate) fn unset_label(&mut self) {
        self.reader.unset_label()
    }

    pub(crate) fn advance(&mut self, amount: usize) {
        self.reader.advance(amount);
        if let Some(bar) = &self.bar {
            bar.inc(amount as u64);
        }
    }

    pub(crate) fn eof(&self) -> bool {
        self.reader.eof()
    }

    pub(crate) fn finish(&mut self) {
        self.reader.finish()
    }

    pub(crate) fn as_slice(&self) -> &'a [u8] {
        self.reader.as_slice()
    }

    pub(crate) fn take_slice(&self, size: usize) -> Option<SliceContent<'a>> {
        self.reader.take_slice(size)
    }
}

impl<'a> SliceReader<'a> {
    pub(crate) fn new(slice: &'a [u8]) -> Self {
        Self {
            slice,
            pos: 0,
            label: None,
        }
    }

    pub(crate) fn label(&self) -> Option<&'static str> {
        self.label
    }

    pub(crate) fn set_label(&mut self, label: &'static str) {
        self.label = Some(label);
    }

    pub(crate) fn unset_label(&mut self) {
        self.label = None;
    }

    pub(crate) fn advance(&mut self, amount: usize) {
        self.pos += amount;
    }

    pub(crate) fn eof(&self) -> bool {
        self.pos >= self.slice.len()
    }

    pub(crate) fn finish(&mut self) {
        self.pos = self.slice.len()
    }

    pub(crate) fn as_slice(&self) -> &'a [u8] {
        &self.slice[self.pos ..]
    }

    pub(crate) fn take_slice(&self, size: usize) -> Option<SliceContent<'a>> {
        // If we've reached or passed the end of the input buffer, stop reading
        if self.eof() {
            return None;
        }
        let chunk;
        let end = self.pos + size;
        // take consideration of the last chunk
        if end >= self.slice.len() {
            // SAFETY: We're within bounds because pos < bytes.len()
            chunk = *unsafe { &self.slice.get_unchecked(self.pos ..) };
            Some(SliceContent::Eof(chunk))
        } else {
            // SAFETY: end < bytes.len() is guaranteed here
            chunk = *unsafe { &self.slice.get_unchecked(self.pos .. end) };
            Some(SliceContent::Data(chunk))
        }
    }
}

pub(crate) enum SliceContent<'a> {
    Data(&'a [u8]), // Normal chunk
    Eof(&'a [u8]),  // Final chunk (possibly partial)
}

impl<'a> SliceContent<'a> {
    pub(crate) fn eof(&self) -> bool {
        match self {
            SliceContent::Data(_) => false,
            SliceContent::Eof(_) => true,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn len(&self) -> usize {
        self.as_slice().len()
    }

    /// Returns the raw bytes of the chunk, if any
    pub(crate) fn as_slice(&self) -> &'a [u8] {
        match self {
            SliceContent::Data(bytes) => bytes,
            SliceContent::Eof(bytes) => bytes,
        }
    }

    pub(crate) fn newlines(&self) -> Vec<usize> {
        match self {
            Self::Data(bytes) => memchr_iter(b'\n', bytes).collect(),
            Self::Eof(bytes) => {
                // Find all newline character positions in the chunk
                let mut newlines: Vec<usize> = memchr_iter(b'\n', bytes).collect();

                // If there are no newlines, treat the whole chunk as a single (newline-less) line
                if newlines.is_empty() {
                    newlines.push(bytes.len());
                } else if (bytes.len() - 1)
                // Ensure the last line is included even if it doesn't end with `\n`
                // SAFETY: newlines is non-empty, so unwrap_unchecked is safe
                != *unsafe { newlines.get_unchecked(newlines.len() - 1) }
                {
                    newlines.push(bytes.len());
                }
                newlines
            }
        }
    }
}

#[derive(Debug)]
pub(crate) struct SliceLineReader<'a, I: Iterator<Item = usize>> {
    reader: SliceReader<'a>,
    breaks: I,
}

impl<'a> SliceLineReader<'a, Memchr<'a>> {
    #[allow(dead_code)]
    pub(crate) fn new(slice: &'a [u8]) -> Self {
        // Pre-compute all newline offsets
        Self::with_breaks(memchr_iter(b'\n', slice).into_iter(), slice)
    }
}

impl<'a, I: Iterator<Item = usize>> SliceLineReader<'a, I> {
    pub(crate) fn with_breaks(breaks: I, slice: &'a [u8]) -> Self {
        Self {
            breaks,
            reader: SliceReader::new(slice),
        }
    }

    pub(crate) fn label(&self) -> Option<&'static str> {
        self.reader.label()
    }

    pub(crate) fn set_label(&mut self, label: &'static str) {
        self.reader.set_label(label)
    }

    #[allow(dead_code)]
    pub(crate) fn unset_label(&mut self) {
        self.reader.unset_label()
    }

    pub(crate) fn read_line(&mut self) -> Option<&'a [u8]> {
        if self.reader.eof() {
            return None;
        }

        if let Some(slice_break) = self.breaks.next() {
            let start = self.reader.pos;
            let end = slice_break;
            self.reader.pos = end + 1; // Move past the newline
            Some(unsafe { &self.reader.slice.get_unchecked(start .. end) })
        } else {
            let bytes = self.reader.as_slice();
            self.reader.finish();
            Some(bytes)
        }
    }

    #[allow(dead_code)]
    pub(crate) fn read_line_inclusive(&mut self) -> Option<&'a [u8]> {
        if self.reader.eof() {
            return None;
        }

        if let Some(slice_break) = self.breaks.next() {
            let start = self.reader.pos;
            let end = slice_break;
            self.reader.pos = end + 1; // Move past the newline
            Some(unsafe { &self.reader.slice.get_unchecked(start ..= end) })
        } else {
            let bytes = self.reader.as_slice();
            self.reader.finish();
            Some(bytes)
        }
    }
}

impl<'a, I: Iterator<Item = usize>> Iterator for SliceLineReader<'a, I> {
    type Item = &'a [u8];
    fn next(&mut self) -> Option<Self::Item> {
        self.read_line()
    }
}

#[cfg(test)]
mod test_reader {
    use super::*;

    // Helper function to create a sample slice reader
    fn create_slice_reader(data: &[u8]) -> SliceReader {
        SliceReader::new(data)
    }

    // Helper function to create a sample SliceLineReader with pre-calculated newlines
    fn create_slice_break_reader(data: &[u8]) -> SliceLineReader<Memchr<'_>> {
        SliceLineReader::new(data)
    }

    #[test]
    fn test_take_slice() {
        let data = b"Hello, world!\nThis is a test.\nAnother line\n";
        let mut reader = create_slice_reader(data);

        // Taking bytes
        if let Some(chunk) = reader.take_slice(10) {
            assert_eq!(chunk.as_slice(), b"Hello, wor");
            reader.advance(chunk.len());
        } else {
            panic!("Expected chunk but got None");
        }

        // Taking another bytes
        if let Some(chunk) = reader.take_slice(10) {
            assert_eq!(chunk.as_slice(), b"ld!\nThis i");
            reader.advance(chunk.len());
        } else {
            panic!("Expected chunk but got None");
        }

        // Taking the final bytes
        if let Some(chunk) = reader.take_slice(100) {
            assert_eq!(chunk.as_slice(), b"s a test.\nAnother line\n");
            reader.advance(chunk.len());
        }

        // No more data, should return None
        assert!(reader.take_slice(100).is_none());
    }

    #[test]
    fn test_slice_break_reader() {
        let data = b"Line 1\nLine 2\nLine 3\nLine 4\n";
        let mut reader = create_slice_break_reader(data);

        // Reading lines
        let line1 = reader.read_line().unwrap();
        assert_eq!(line1, b"Line 1");

        let line2 = reader.read_line().unwrap();
        assert_eq!(line2, b"Line 2");

        let line3 = reader.read_line().unwrap();
        assert_eq!(line3, b"Line 3");

        let line4 = reader.read_line().unwrap();
        assert_eq!(line4, b"Line 4");

        // Should return None when EOF is reached
        assert!(reader.read_line().is_none());
    }

    #[test]
    fn test_slice_break_reader_inclusive() {
        let data = b"Line 1\nLine 2\nLine 3\n";
        let mut reader = create_slice_break_reader(data);

        // Reading lines inclusively
        let line1 = reader.read_line_inclusive().unwrap();
        assert_eq!(line1, b"Line 1\n");

        let line2 = reader.read_line_inclusive().unwrap();
        assert_eq!(line2, b"Line 2\n");

        let line3 = reader.read_line_inclusive().unwrap();
        assert_eq!(line3, b"Line 3\n");

        // Should return None when EOF is reached
        assert!(reader.read_line_inclusive().is_none());
    }

    #[test]
    fn test_progress_bar_integration() {
        let data = b"Hello\nWorld\nTest\n";

        let pb = ProgressBar::new(3);
        let mut reader_with_bar = SliceProgressBarReader::new(data);

        // Simulate chunk reading with progress bar updates
        while let Some(chunk) = reader_with_bar.take_slice(10) {
            reader_with_bar.advance(chunk.len());
        }
        pb.finish_with_message("Done reading chunks");
    }

    #[test]
    fn test_no_label_assigned() {
        let data = b"Test data with no label";
        let reader = create_slice_reader(data);

        assert!(reader.label().is_none());
    }

    #[test]
    fn test_label_assigned() {
        let data = b"Test data with label";
        let mut reader = create_slice_reader(data);

        reader.set_label("MyLabel");
        assert_eq!(reader.label(), Some("MyLabel"));
    }

    #[test]
    fn test_take_slice_empty() {
        let data = b"";
        let reader = create_slice_reader(data);

        // Should return None for an empty slice
        assert!(reader.take_slice(10).is_none());
    }
}
