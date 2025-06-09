use anyhow::{Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};

pub struct ChunkIO<IO, I, O>
where
    IO: ChunkIOFactory<I, O>,
    I: ChunkInput,
    O: ChunkOutput<Chunk = I::Chunk>,
{
    io: IO,
    n_queue: usize,
    _i: std::marker::PhantomData<I>,
    _o: std::marker::PhantomData<O>,
}

impl<IO, I, O> ChunkIO<IO, I, O>
where
    IO: ChunkIOFactory<I, O>,
    I: ChunkInput,
    O: ChunkOutput<Chunk = I::Chunk>,
{
    pub fn new(io: IO, n_queue: usize) -> Self {
        Self {
            io,
            n_queue,
            _i: std::marker::PhantomData,
            _o: std::marker::PhantomData,
        }
    }

    /// Spawn a thread to run the **output logic** (output must be `Send`)
    pub fn chunk_io_spawn_output(&self) -> Result<O::Output>
    where
        O::Output: Send,
    {
        let (tx, rx) = bounded(self.n_queue);
        std::thread::scope(|scope| {
            let io = &self.io;
            // Spawn a new thread to run the output consumer
            // Closures can capture values from their environment in three ways,
            // which directly map to the three ways a function can take a parameter:
            // borrowing immutably, borrowing mutably, and taking ownership. The
            // closure will decide which of these to use based on what the body of
            // the function does with the captured values.
            let handle = scope.spawn(move || -> Result<O::Output> {
                io.new_output().chunk_receive(rx).context("Output error")
            });

            // Run input on the main thread
            io.new_input().chunk_send(tx).context("Input error")?;

            // Wait for the output thread to finish and return its result
            match handle.join() {
                Ok(inner) => inner,
                Err(e) => Err(anyhow::anyhow!(
                    "Thread panicked during output processing: {:?}",
                    e
                )),
            }
        })
    }

    /// Spawn a thread to run the **input logic**
    pub fn chunk_io_spawn_input(&self) -> Result<O::Output> {
        let (tx, rx) = bounded(self.n_queue);
        std::thread::scope(|scope| -> Result<O::Output> {
            let io = &self.io;
            // Spawn a new thread to run the input producer
            // we use move here to ensure the tx will be dropped when closed
            let handle = scope.spawn(move || -> Result<()> {
                io.new_input()
                    .chunk_send(tx)
                    .context("`chunk_send()` error")
            });

            // Run output on the main thread
            let out = io
                .new_output()
                .chunk_receive(rx)
                .context("`chunk_receive()` error")?;

            // Wait for input thread and propagate errors if any
            handle.join().map_err(|e| {
                anyhow::anyhow!(
                    "Thread panicked during input processing: {:?}",
                    e
                )
            })?;
            Ok(out)
        })
    }
}

pub trait ChunkIOFactory<I, O>: Sync
where
    I: ChunkInput,
    O: ChunkOutput<Chunk = I::Chunk>,
{
    fn new_input(&self) -> I;
    fn new_output(&self) -> O;
}

pub trait ChunkInput {
    type Chunk: Send;

    fn chunk_send(self, tx: Sender<Self::Chunk>) -> Result<()>;
}

pub trait ChunkOutput {
    type Output;
    type Chunk: Send;
    // Do something for the chunk
    fn chunk_receive(self, rx: Receiver<Self::Chunk>) -> Result<Self::Output>;
}
