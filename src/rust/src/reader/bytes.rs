use std::{
    cell::{Cell, RefCell},
    io::Read,
};

use bytes::{Bytes, BytesMut};
use indicatif::ProgressBar;
use memchr::{memchr, memchr_iter};

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

#[allow(dead_code)]
#[derive(Debug)]
pub(crate) struct BytesReader<R: Read> {
    reader: R,
    leftover: BytesMut,
    label: Option<&'static str>,
}

impl<R: Read> BytesReader<R> {
    pub(crate) fn new(reader: R) -> Self {
        Self {
            reader,
            leftover: BytesMut::new(),
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

    pub(crate) fn borrow_reader(&mut self) -> &mut R {
        self.reader.by_ref()
    }

    pub(crate) fn read_bytes(&mut self, size: usize) -> std::io::Result<Option<BytesContent>> {
        let leftover_len = self.leftover.len();
        let mut buf;
        if self.leftover.len() > 0 {
            buf = std::mem::take(&mut self.leftover);
            buf.reserve(size);
        } else {
            buf = BytesMut::with_capacity(size);
        }
        unsafe { buf.set_len(leftover_len + size) };
        let nbytes = self.reader.read(&mut buf[leftover_len ..])?;
        unsafe { buf.set_len(leftover_len + nbytes) };

        if buf.len() == 0 {
            Ok(None)
        } else if nbytes == 0 {
            Ok(Some(BytesContent::Eof(buf)))
        } else {
            Ok(Some(BytesContent::Data(buf)))
        }
    }

    pub(crate) fn take_leftover(&mut self, leftover: BytesMut) {
        self.leftover = leftover;
    }
}

pub(crate) enum BytesContent {
    Data(BytesMut), // Normal chunk
    Eof(BytesMut),  // Final chunk (possibly partial)
}

impl BytesContent {
    pub(crate) fn eof(&self) -> bool {
        match self {
            BytesContent::Data(_) => false,
            BytesContent::Eof(_) => true,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn len(&self) -> usize {
        self.borrow_bytes().len()
    }

    #[allow(dead_code)]
    pub(crate) fn borrow_bytes(&self) -> &BytesMut {
        match self {
            BytesContent::Data(bytes) => bytes,
            BytesContent::Eof(bytes) => bytes,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn borrow_mut_bytes(&mut self) -> &mut BytesMut {
        match self {
            BytesContent::Data(ref mut bytes) => bytes,
            BytesContent::Eof(ref mut bytes) => bytes,
        }
    }

    pub(crate) fn into_bytes(self) -> BytesMut {
        match self {
            BytesContent::Data(bytes) => bytes,
            BytesContent::Eof(bytes) => bytes,
        }
    }

    pub(crate) fn as_slice(&self) -> &[u8] {
        &self.borrow_bytes()
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

// Used by Fastq reader, it will return borrowed line slice instead
#[derive(Debug)]
pub(crate) struct BytesBreaksReader<I: Iterator<Item = usize>> {
    bytes: Bytes,
    breaks: RefCell<I>,
    offset: Cell<usize>,
    label: Option<&'static str>,
}

impl<'a> BytesBreaksReader<std::vec::IntoIter<usize>> {
    #[allow(dead_code)]
    pub(crate) fn new(bytes: Bytes) -> Self {
        // Pre-compute all newline offsets
        let breaks = memchr_iter(b'\n', &bytes)
            .into_iter()
            .collect::<Vec<usize>>()
            .into_iter();
        Self::with_breaks(breaks, bytes)
    }
}

impl<I: Iterator<Item = usize>> BytesBreaksReader<I> {
    pub(crate) fn with_breaks(breaks: I, bytes: Bytes) -> Self {
        Self {
            breaks: RefCell::new(breaks),
            bytes,
            offset: Cell::new(0),
            label: None,
        }
    }

    pub(crate) fn label(&self) -> Option<&'static str> {
        self.label
    }

    pub(crate) fn set_label(&mut self, label: &'static str) {
        self.label = Some(label);
    }

    #[allow(dead_code)]
    pub(crate) fn unset_label(&mut self) {
        self.label = None
    }

    pub(crate) fn borrow_bytes(&self) -> &Bytes {
        &self.bytes
    }

    #[allow(dead_code)]
    pub(crate) fn as_slice(&self) -> &[u8] {
        &self.borrow_bytes()
    }

    pub(crate) fn read_line(&self) -> Option<&[u8]> {
        let offset = self.offset.get();
        if offset >= self.bytes.len() {
            return None;
        }

        if let Some(bytes_break) = self.breaks.borrow_mut().next() {
            let end = bytes_break;
            self.offset.set(end + 1); // Move past the newline
            Some(unsafe { &self.bytes.get_unchecked(offset .. end) })
        } else {
            let bytes = &self.bytes[offset ..];
            self.offset.set(bytes.len());
            Some(bytes)
        }
    }

    #[allow(dead_code)]
    pub(crate) fn read_line_inclusive(&self) -> Option<&[u8]> {
        let offset = self.offset.get();
        if offset >= self.bytes.len() {
            return None;
        }

        if let Some(bytes_break) = self.breaks.borrow_mut().next() {
            let end = bytes_break;
            self.offset.set(end + 1); // Move past the newline
            Some(unsafe { &self.bytes.get_unchecked(offset ..= end) })
        } else {
            let bytes = &self.bytes[offset ..];
            self.offset.set(bytes.len());
            Some(bytes)
        }
    }
}

#[derive(Debug)]
pub(crate) struct BytesLineReader {
    bytes: Bytes,
    offset: usize,
    label: Option<&'static str>,
}

impl BytesLineReader {
    pub(crate) fn new(bytes: Bytes) -> Self {
        Self {
            bytes,
            offset: 0,
            label: None,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn label(&self) -> Option<&'static str> {
        self.label
    }

    #[allow(dead_code)]
    pub(crate) fn set_label(&mut self, label: &'static str) {
        self.label = Some(label);
    }

    #[allow(dead_code)]
    pub(crate) fn unset_label(&mut self) {
        self.label = None
    }

    #[allow(dead_code)]
    pub(crate) fn borrow_bytes(&self) -> &Bytes {
        &self.bytes
    }

    #[allow(dead_code)]
    pub(crate) fn as_slice(&self) -> &[u8] {
        &self.borrow_bytes()
    }

    #[allow(dead_code)]
    pub(crate) fn read_line(&mut self) -> Option<Bytes> {
        if self.offset >= self.bytes.len() {
            return None;
        }
        let line;
        if let Some(pos) = memchr(b'\n', &self.bytes[self.offset ..]) {
            line = unsafe { self.bytes.get_unchecked(self.offset .. self.offset + pos) };
            self.offset += pos + 1; // Move past the newline
        } else {
            line = unsafe { self.bytes.get_unchecked(self.offset ..) };
            self.offset = self.bytes.len();
        }
        Some(self.bytes.slice_ref(line))
    }

    pub(crate) fn read_line_inclusive(&mut self) -> Option<Bytes> {
        if self.offset >= self.bytes.len() {
            return None;
        }
        let line;
        if let Some(pos) = memchr(b'\n', &self.bytes[self.offset ..]) {
            line = unsafe { self.bytes.get_unchecked(self.offset ..= self.offset + pos) };
            self.offset += pos + 1; // Move past the newline
        } else {
            line = unsafe { self.bytes.get_unchecked(self.offset ..) };
            self.offset = self.bytes.len();
        }
        Some(self.bytes.slice_ref(line))
    }
}
