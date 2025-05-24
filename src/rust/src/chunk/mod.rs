use std::fmt::Display;
use std::fs::File;
use std::io::BufWriter;

mod io;
mod reader;
mod splitter;

use io::ChunkIO;
use reader::ChunkReader;
use splitter::ChunkSplitter;

struct ChunkFactory<F>
where
    F: ChunkProcessor,
{
    processor: F,
    buffersize: usize,
}

impl<F> ChunkFactory<F>
where
    F: ChunkProcessor,
{
    pub fn new(processor: F, buffersize: usize) -> Self {
        Self {
            processor,
            buffersize,
        }
    }
}

pub trait ChunkParser {
    fn parse<F>(&self, chunk: Vec<u8>, push: &mut F) -> Result<(), String>
    where
        F: FnMut(Vec<u8>) -> Result<(), String>;
}

pub trait ChunkProcessor: Sync + Sized {
    type Parser: ChunkParser;

    fn new_parser(&self) -> Self::Parser;

    fn new_splitter(&self) -> ChunkSplitter {
        ChunkSplitter::Byte(b'\n')
    }

    #[allow(clippy::too_many_arguments)]
    fn chunk_io<P>(
        self,
        input: P,
        output: P,
        read_buffer: usize,
        write_buffer: usize,
        parse_buffer: usize,
        read_queue: usize,
        write_queue: usize,
        threads: usize,
    ) -> std::result::Result<(), String>
    where
        P: AsRef<std::path::Path> + Display,
    {
        // one thread is kept for writer
        let reader = ChunkReader::build(
            File::open(&input).map_err(|e| {
                format!("Cannot open {}: {}", input.as_ref().display(), e)
            })?,
            self.new_splitter(),
            read_buffer,
        );
        let writer = BufWriter::with_capacity(
            write_buffer,
            File::create(&output).map_err(|e| {
                format!(
                    "Failed to create file {}: {}",
                    output.as_ref().display(),
                    e
                )
            })?,
        );
        let factory = ChunkFactory::new(self, parse_buffer);
        let mut io = ChunkIO::new(
            reader,
            writer,
            factory,
            read_queue,
            write_queue,
            threads,
        );
        io.process()
    }
}
