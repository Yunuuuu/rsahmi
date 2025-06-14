use std::io::{self, Read};
use std::mem::MaybeUninit;

pub struct Buffer {
    // The buffer (partially initialized).
    buf: Box<[MaybeUninit<u8>]>,
    // The current seek offset into `buf`, must always be <= `filled`.
    pos: usize,
    // How many bytes are currently filled (valid data).
    filled: usize,
    // How many bytes are initialized (used for safety).
    initialized: usize,
}

impl Buffer {
    #[inline]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            buf: Box::new_uninit_slice(capacity),
            pos: 0,
            filled: 0,
            initialized: 0,
        }
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        self.buf.len()
    }

    #[inline]
    #[allow(dead_code)]
    pub fn pos(&self) -> usize {
        self.pos
    }

    #[inline]
    #[allow(dead_code)]
    pub fn filled(&self) -> usize {
        self.filled
    }

    #[inline]
    pub fn buffer(&self) -> &[u8] {
        debug_assert!(self.pos <= self.filled);
        unsafe {
            std::slice::from_raw_parts(
                self.buf.as_ptr().add(self.pos) as *const u8,
                self.filled - self.pos,
            )
        }
    }

    #[inline]
    pub fn extend(&mut self, offset: usize) {
        if offset == 0 {
            return;
        }
        let mut new_buf = Box::new_uninit_slice(self.buf.len() + offset);

        // SAFETY: Only copying the filled part
        unsafe {
            std::ptr::copy_nonoverlapping(
                self.buf.as_ptr(),
                new_buf.as_mut_ptr(),
                self.filled,
            );
        }

        self.buf = new_buf;
        self.initialized = self.filled;
    }

    #[inline]
    pub fn consume(&mut self, amt: usize) {
        self.pos = std::cmp::min(self.pos + amt, self.filled);
    }

    #[inline]
    #[allow(dead_code)]
    pub fn discard(&mut self) {
        self.pos = 0;
        self.filled = 0;
    }

    // Used internally to remove the consumed buffer
    #[inline]
    pub fn compact(&mut self) {
        // Compact the buffer if there is consumed space at the beginning.
        debug_assert!(self.pos <= self.filled);
        if self.pos > 0 {
            if self.pos == self.filled {
                self.discard();
            } else {
                // SAFETY: the region [pos..filled] is initialized and valid
                // Thought the two ptr are overlapped, we don't need
                // use copy directly, since we won't use the source bytes after copying
                unsafe {
                    std::ptr::copy(
                        self.buf.as_ptr().add(self.pos),
                        self.buf.as_mut_ptr(),
                        self.filled - self.pos,
                    );
                }
                self.filled = self.filled - self.pos;
                self.pos = 0;
            }
        }
    }

    /// Fill the buffer
    pub fn fill_buf(
        &mut self,
        mut reader: impl Read,
    ) -> io::Result<Option<usize>> {
        let capacity = self.capacity();
        if self.filled >= capacity {
            return Ok(Some(0)); // No space left
        }

        let remaining = capacity - self.filled;
        // SAFETY: We're only constructing a mutable slice to the un-filled portion.
        let unfilled_slice = unsafe {
            std::slice::from_raw_parts_mut(
                self.buf.as_mut_ptr().add(self.filled) as *mut u8,
                remaining,
            )
        };

        let n = reader.read(unfilled_slice)?;
        if n == 0 {
            return Ok(None);
        }

        self.filled += n;
        self.initialized = self.initialized.max(self.filled);

        Ok(Some(n))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_buffer_initialization() {
        let buf = Buffer::with_capacity(128);
        assert_eq!(buf.capacity(), 128);
        assert_eq!(buf.pos(), 0);
        assert_eq!(buf.filled(), 0);
    }

    #[test]
    fn test_fill_buf_and_borrow() {
        let mut buffer = Buffer::with_capacity(16);
        let data = b"hello world";
        let mut cursor = Cursor::new(&data[..]);

        let n = buffer.fill_buf(&mut cursor).unwrap().unwrap();
        assert_eq!(n, data.len());

        let borrowed = buffer.buffer();
        assert_eq!(borrowed, b"hello world");
    }

    #[test]
    fn test_consume_and_compact() {
        let mut buffer = Buffer::with_capacity(16);
        let data = b"0123456789";
        let mut cursor = Cursor::new(&data[..]);

        buffer.fill_buf(&mut cursor).unwrap().unwrap();
        assert_eq!(buffer.buffer(), b"0123456789");

        buffer.consume(4);
        assert_eq!(buffer.buffer(), b"456789");

        buffer.compact();
        assert_eq!(buffer.pos(), 0);
        assert_eq!(buffer.buffer(), b"456789");

        // ensure buffer can still fill more data
        let data2 = b"ABCD";
        let mut cursor2 = Cursor::new(&data2[..]);
        buffer.fill_buf(&mut cursor2).unwrap().unwrap();
        let result = buffer.buffer();
        assert!(result.ends_with(b"ABCD"));
    }

    #[test]
    fn test_extend() {
        let mut buffer = Buffer::with_capacity(8);
        let data = b"abc";
        let mut cursor = Cursor::new(&data[..]);

        buffer.fill_buf(&mut cursor).unwrap().unwrap();
        assert_eq!(buffer.buffer(), b"abc");

        buffer.extend(8);
        assert!(buffer.capacity() >= 16); // new buffer is larger
        assert_eq!(buffer.buffer(), b"abc"); // old data retained
    }

    #[test]
    fn test_fill_buf_eof() {
        let mut buffer = Buffer::with_capacity(16);
        let empty = Cursor::new(b"");
        let r = buffer.fill_buf(empty).unwrap();
        assert_eq!(r, None);
        assert_eq!(buffer.buffer().len(), 0);
    }

    #[test]
    fn test_discard() {
        let mut buffer = Buffer::with_capacity(16);
        let mut cursor = Cursor::new(b"test data");
        buffer.fill_buf(&mut cursor).unwrap().unwrap();
        buffer.consume(4);
        assert!(buffer.buffer().starts_with(b" data"));

        buffer.discard();
        assert_eq!(buffer.filled(), 0);
        assert_eq!(buffer.pos(), 0);
        assert_eq!(buffer.buffer().len(), 0);
    }
}
