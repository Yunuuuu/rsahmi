use anyhow::{Context, Result};
use crossbeam_channel::{bounded, Receiver, Sender};

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

/// Spawn a thread to run the **output logic** (output must be `Send`)
pub fn chunk_io_spawn_output<I, O>(
    input: I,
    output: O,
    n_queue: usize,
) -> Result<()>
where
    I: ChunkInput,
    O: ChunkOutput<Chunk = I::Chunk> + Send,
{
    let (tx, rx) = bounded::<I::Chunk>(n_queue);
    std::thread::scope(|scope| {
        // Spawn a new thread to run the output consumer
        // Closures can capture values from their environment in three ways,
        // which directly map to the three ways a function can take a parameter:
        // borrowing immutably, borrowing mutably, and taking ownership. The
        // closure will decide which of these to use based on what the body of
        // the function does with the captured values.
        let handle = scope.spawn(move || -> Result<()> {
            output.chunk_receive(rx).context("Output error")
        });

        // Run input on the main thread
        input.chunk_send(tx).context("Input error")?;

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

/// Spawn a thread to run the **input logic** (input must be `Send`)
pub fn chunk_io_spawn_input<I, O>(
    input: I,
    output: O,
    n_queue: usize,
) -> Result<()>
where
    I: ChunkInput + Send,
    O: ChunkOutput<Chunk = I::Chunk>,
{
    let (tx, rx) = bounded::<I::Chunk>(n_queue);
    std::thread::scope(|scope| {
        // Spawn a new thread to run the input producer
        // we use move here to ensure the tx will be dropped when closed
        let handle = scope.spawn(move || -> Result<()> {
            input.chunk_send(tx).context("`chunk_send()` error")
        });

        // Run output on the main thread
        output
            .chunk_receive(rx)
            .context("`chunk_receive()` error")?;

        // Wait for input thread and propagate errors if any
        match handle.join() {
            Ok(inner) => inner,
            Err(e) => Err(anyhow::anyhow!(
                "Thread panicked during input processing: {:?}",
                e
            )),
        }
    })
}
#[cfg(test)]
mod tests {
    use std::cell::RefCell;

    use super::*;

    struct DummyInput {
        data: Vec<i32>,
    }

    impl ChunkInput for DummyInput {
        type Chunk = i32;

        fn chunk_send(self, tx: Sender<Self::Chunk>) -> Result<()> {
            for val in self.data {
                tx.send(val)?;
            }
            Ok(())
        }
    }

    struct DummyOutput {
        collected: Vec<i32>,
    }

    impl ChunkOutput for DummyOutput {
        type Chunk = i32;

        fn chunk_receive(mut self, rx: Receiver<Self::Chunk>) -> Result<()> {
            for chunk in rx.iter() {
                self.collected.push(chunk);
            }
            Ok(())
        }
    }

    #[test]
    fn test_chunk_io_spawn_output() -> Result<()> {
        let input = DummyInput {
            data: vec![1, 2, 3, 4],
        };
        let output = DummyOutput {
            collected: Vec::new(),
        };

        chunk_io_spawn_output(input, output, 2)?;

        // `output` was moved, so we can't check it here.
        // Instead, test output inside `chunk_receive` or use an external collector if needed.
        Ok(())
    }

    #[test]
    fn test_chunk_io_spawn_input() -> Result<()> {
        let input = DummyInput {
            data: vec![10, 20, 30],
        };
        let mut collected = Vec::new();

        let output = DummyOutput {
            collected: RefCell::new(&mut collected),
        };

        chunk_io_spawn_input(input, output, 2)?;
        assert_eq!(collected, vec![10, 20, 30]);

        Ok(())
    }
}
