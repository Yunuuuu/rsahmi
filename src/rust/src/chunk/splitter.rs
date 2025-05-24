pub enum ChunkSplitter {
    Byte(u8),           // Split at a specific byte, like `b'\n'`
    FixedLength(usize), // Split every N bytes
}

impl ChunkSplitter {
    pub fn breakpoint(&self, buffer: &[u8]) -> Option<usize> {
        match self {
            ChunkSplitter::Byte(b) => memchr::memrchr(*b, buffer),
            ChunkSplitter::FixedLength(n) => {
                if buffer.len() >= *n {
                    Some(*n - 1)
                } else {
                    None
                }
            }
        }
    }
}

impl Default for ChunkSplitter {
    fn default() -> Self {
        ChunkSplitter::Byte(b'\n') // Default to splitting at newline
    }
}
