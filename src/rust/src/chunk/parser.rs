use anyhow::{Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};

use crate::chunk::io::{ChunkInput, ChunkOutput};

pub struct ChunkParallelIO<IO, F, P>
where
    F: ChunkParserFactory<Parser = P>,
    P: ChunkParser,
{
    io: IO,
    // sender: Sender<F::Chunk>,
    // receiver: Receiver<F::Chunk>,
    factory: F,
    threads: usize,
    n_queue: usize,
}

pub trait ChunkParserFactory: Sync {
    type Parser: ChunkParser;
    fn new_parser(&self) -> Self::Parser;
}

pub trait ChunkParser {
    type Chunk: Send;
    type Record: Send;
    fn parse(
        &self,
        rx: Receiver<Self::Chunk>,
        tx: Sender<Self::Record>,
    ) -> Result<()>;
}

impl<IO, F, P> ChunkParallelIO<IO, F, P>
where
    F: ChunkParserFactory<Parser = P>,
    P: ChunkParser,
{
    pub fn from_input(
        io: IO,
        factory: F,
        threads: usize,
        n_queue: usize,
    ) -> Self
    where
        IO: ChunkInput,
    {
        Self {
            io,
            factory,
            threads,
            n_queue,
        }
    }
    pub fn from_output(
        io: IO,
        factory: F,
        threads: usize,
        n_queue: usize,
    ) -> Self
    where
        IO: ChunkOutput,
    {
        Self {
            io,
            factory,
            threads,
            n_queue,
        }
    }
}

// for Input decorator, we send the parsed record
impl<I, F, P> ChunkInput for ChunkParallelIO<I, F, P>
where
    I: ChunkInput<Chunk = P::Chunk> + Send,
    F: ChunkParserFactory<Parser = P>,
    P: ChunkParser,
{
    type Chunk = P::Record;
    fn chunk_send(self, tx: Sender<Self::Chunk>) -> Result<()> {
        let (parser_tx, parser_rx) = bounded::<P::Chunk>(self.n_queue);
        std::thread::scope(|scope| {
            // Spawn input thread
            let io = self.io;
            let handle =
                scope.spawn(move || -> Result<()> { io.chunk_send(parser_tx) });

            // Spawn parser threads
            let mut handle_list = Vec::with_capacity(self.threads);
            for _ in 0 .. self.threads {
                let rx = parser_rx.clone();
                let tx = tx.clone();
                let factory = &self.factory;
                let handle = scope.spawn(move || {
                    let parser = factory.new_parser();
                    parser.parse(rx, tx)
                });
                handle_list.push(handle);
            }
            handle
                .join()
                .map_err(|e| anyhow::anyhow!("{:?}", e))
                .context("Input thread panicked")?
                .context("Inner `chunk_send()` failed")?;
            for handle in handle_list {
                handle
                    .join()
                    .map_err(|e| anyhow::anyhow!("{:?}", e))
                    .context("Parser thread panicked")?
                    .context("`parse()` failed")?;
            }
            Ok(())
        })
    }
}

// for Output decorator, we send the parsed record
impl<O, F, P> ChunkOutput for ChunkParallelIO<O, F, P>
where
    O: ChunkOutput<Chunk = P::Record> + Send,
    F: ChunkParserFactory<Parser = P>,
    P: ChunkParser,
{
    type Chunk = P::Record;
    fn chunk_receive(self, rx: Receiver<Self::Chunk>) -> Result<()> {}
}
