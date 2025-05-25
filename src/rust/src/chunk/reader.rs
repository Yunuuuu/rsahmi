use std::fs::File;
use std::io::Read;
use std::path::Path;

use super::splitter::ChunkSplitter;

pub struct ChunkReader<R>
where
    R: Read,
{
    // reader to read from
    reader: R,
    // size of the buffer to read
    buffersize: usize,
    // buffer to store lines
    buffer: Vec<u8>,
    // used to store the byte from where to read
    offset: usize,
    splitter: ChunkSplitter,
}

impl<R> Default for ChunkReader<R>
where
    R: Read + Default,
{
    fn default() -> Self {
        Self::new(R::default())
    }
}

impl<P> From<P> for ChunkReader<File>
where
    P: AsRef<Path>,
{
    fn from(path: P) -> Self {
        let reader = File::open(&path).unwrap_or_else(|e| {
            panic!("Failed to open file {}: {}", path.as_ref().display(), e)
        });
        Self::new(reader)
    }
}

impl<R> ChunkReader<R>
where
    R: Read,
{
    pub fn new(reader: R) -> Self {
        Self::build(reader, ChunkSplitter::default(), 1024 * 1024)
    }

    pub fn build(
        reader: R,
        splitter: ChunkSplitter,
        buffersize: usize,
    ) -> Self {
        Self {
            reader,
            buffersize,
            buffer: vec![0; buffersize],
            offset: 0,
            splitter,
        }
    }

    fn take_next(&mut self) -> Option<std::result::Result<Vec<u8>, String>> {
        match self.read_buffer() {
            Ok(0) => self.take_leftover().map(Ok),
            Ok(_) => match self.take_chunk() {
                Some(chunk) => Some(Ok(chunk)),
                None => {
                    // if we have not found a chunk, we need to read more
                    if self.is_buffer_full() {
                        if let Err(e) = self.extend_buffer() {
                            return Some(Err(e));
                        }
                    }
                    self.take_next()
                }
            },
            Err(e) => Some(Err(e)),
        }
    }

    // read data from the reader into the buffer
    fn read_buffer(&mut self) -> std::result::Result<usize, String> {
        let read_bytes = self
            .reader
            .read(&mut self.buffer[self.offset ..])
            .map_err(|e| format!("Read failed: {e}"))?;
        self.offset += read_bytes;
        Ok(read_bytes)
    }

    fn take_chunk(&mut self) -> Option<Vec<u8>> {
        self.splitter
            .breakpoint(&self.buffer[.. self.offset])
            .map(|pos| self.take_buffer(pos))
    }

    fn take_leftover(&mut self) -> Option<Vec<u8>> {
        // ensure all buffer has been read
        if self.offset > 0 {
            let mut chunk = self.take_buffer(self.offset - 1);
            // always ensure the last line has `\n`
            if !chunk.ends_with(b"\n") {
                chunk.push(b'\n');
            }
            Some(chunk)
        } else {
            None
        }
    }

    fn take_buffer(&mut self, pos: usize) -> Vec<u8> {
        let chunk = self.buffer[..= pos].to_vec(); // copy out the chunk
        let remaining = self.offset - (pos + 1); // move remaining bytes to the front
        if remaining > 0 {
            self.buffer.copy_within((pos + 1) .. self.offset, 0);
        }
        self.offset = remaining;
        chunk
    }

    fn extend_buffer(&mut self) -> std::result::Result<(), String> {
        self.buffer
            .try_reserve(self.buffersize)
            .map_err(|e| format!("Increase buffer failed: {e}"))?;
        // everytime we increase the buffer, we need to fill it with 0
        self.buffer.resize(self.buffer.capacity(), 0);
        Ok(())
    }

    fn is_buffer_full(&self) -> bool {
        self.offset == self.buffer.capacity()
    }
}

impl<R> Iterator for ChunkReader<R>
where
    R: Read,
{
    type Item = std::result::Result<Vec<u8>, String>;

    fn next(&mut self) -> Option<Self::Item> {
        self.take_next()
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_chunk_reader_with_newlines() {
        let data = b"line1\nline2\nline3\n";
        let cursor = Cursor::new(data.to_vec());

        let mut reader = ChunkReader::new(cursor);
        let chunks = reader.next().unwrap().unwrap();

        assert_eq!(chunks.len(), 18);
        assert_eq!(chunks, data);
    }

    #[test]
    fn test_chunk_reader_with_incomplete_final_line() {
        let data = b"line1\nline2\nincomplete";
        let cursor = Cursor::new(data.to_vec());

        let reader = ChunkReader::new(cursor);
        let chunks: Vec<_> = reader.map(Result::unwrap).collect();
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0], b"line1\nline2\n");
        assert_eq!(chunks[1], b"incomplete\n");
    }

    #[test]
    fn test_chunk_reader_empty_input() {
        let data = b"";
        let cursor = Cursor::new(data.to_vec());

        let mut reader = ChunkReader::new(cursor);
        assert!(reader.next().is_none());
    }
}
