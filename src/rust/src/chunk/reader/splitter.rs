use memchr::memmem;

#[allow(dead_code)]
pub enum ChunkSplitter {
    Byte(u8),             // Split at a specific byte, like `b'\n'`
    Bytes(&'static [u8]), // Special case for FASTQ format
}

impl ChunkSplitter {
    pub fn breakpoint(&self, buffer: &[u8]) -> Option<usize> {
        match self {
            ChunkSplitter::Byte(b) => memchr::memrchr(*b, buffer),
            ChunkSplitter::Bytes(bytes) => {
                // FASTQ format has a specific structure, we look for the '@' character
                // which indicates the start of a new record.
                memmem::FinderRev::new(bytes).rfind(buffer)
            }
        }
    }
}

impl Default for ChunkSplitter {
    fn default() -> Self {
        ChunkSplitter::Byte(b'\n') // Default to splitting at newline
    }
}
