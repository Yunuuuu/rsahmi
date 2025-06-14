use anyhow::{anyhow, Result};
use crossbeam_channel::{Receiver, Sender};

// use super::reader::ChunkReader;
// use super::ChunkFactory;
// use super::{ChunkParser, ChunkProcessor};

pub trait ChunkReader {
    type Output: Send;
    fn read(self) -> Result<()>;
}

pub trait ChunkParser {
    type Input: Send;
    type Output: Send;
    fn parse(self) -> Result<()>;
}

pub trait ChunkWriter {
    type Input;
    fn write(self) -> Result<()>;
}

/// Factory trait: responsible for constructing a reader, parser, and writer.
pub trait ChunkFactory<'a>: Sync {
    type Reader: ChunkReader;
    type Parser: ChunkParser<Input = <Self::Reader as ChunkReader>::Output>;
    type Writer: ChunkWriter<Input = <Self::Parser as ChunkParser>::Output>;

    fn new_reader(
        &'a self,
        tx: Sender<<Self::Reader as ChunkReader>::Output>,
    ) -> Result<Self::Reader>;

    fn new_parser(
        &'a self,
        rx: Receiver<<Self::Parser as ChunkParser>::Input>,
        tx: Sender<<Self::Parser as ChunkParser>::Output>,
    ) -> Result<Self::Parser>;

    fn new_writer(
        &'a self,
        rx: Receiver<<Self::Writer as ChunkWriter>::Input>,
    ) -> Result<Self::Writer>;

    #[allow(clippy::type_complexity)]
    fn channel_reader_parser(
        &'a self,
    ) -> (
        Sender<<Self::Reader as ChunkReader>::Output>,
        Receiver<<Self::Reader as ChunkReader>::Output>,
    );

    #[allow(clippy::type_complexity)]
    fn channel_parser_writer(
        &'a self,
    ) -> (
        Sender<<Self::Parser as ChunkParser>::Output>,
        Receiver<<Self::Parser as ChunkParser>::Output>,
    );
}
pub struct ChunkIO<'a, IO>
where
    IO: ChunkFactory<'a>,
{
    io: IO,
    _phantom: std::marker::PhantomData<&'a ()>,
}

impl<'a, IO> ChunkIO<'a, IO>
where
    IO: ChunkFactory<'a>,
{
    pub fn new(io: IO) -> Self {
        Self {
            io,
            _phantom: std::marker::PhantomData,
        }
    }

    pub fn run(&'a self) -> Result<()> {
        let (reader_tx, parser_rx) = self.io.channel_reader_parser();
        let (parser_tx, writer_rx) = self.io.channel_parser_writer();

        // Start the processor threads
        std::thread::scope(move |scope| {
            // read files and pass chunks to parser
            let reader_handle =
                scope.spawn(move || self.io.new_reader(reader_tx)?.read());

            // parser threads, will send each passed lines to `writer`
            let parser_handle = scope.spawn(move || {
                self.io.new_parser(parser_rx, parser_tx)?.parse()
            });

            // Writer thread: write data to the file
            // accept lines from `parser_tx`
            let writer_handle =
                scope.spawn(move || self.io.new_writer(writer_rx)?.write());

            // Collect results from each thread
            reader_handle
                .join()
                .map_err(|e| anyhow!("Reader thread panicked: {:?}", e))??;
            parser_handle
                .join()
                .map_err(|e| anyhow!("Parser thread panicked: {:?}", e))??;
            writer_handle
                .join()
                .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))??;
            Ok(())
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::{Arc, Mutex};

    use anyhow::Result;
    use crossbeam_channel::{bounded, Receiver, Sender};

    use super::*;

    #[derive(Clone)]
    struct MockIO {
        // Shared state to verify the final output
        pub collected: Arc<Mutex<Vec<String>>>,
    }

    struct MockReader {
        tx: Sender<String>,
    }

    impl ChunkReader for MockReader {
        type Output = String;

        fn read(self) -> Result<()> {
            self.tx.send("line1".to_string())?;
            self.tx.send("line2".to_string())?;
            Ok(())
        }
    }

    struct MockParser {
        rx: Receiver<String>,
        tx: Sender<String>,
    }

    impl ChunkParser for MockParser {
        type Input = String;
        type Output = String;

        fn parse(self) -> Result<()> {
            for line in self.rx.iter() {
                let parsed = format!("parsed:{}", line);
                self.tx.send(parsed)?;
            }
            Ok(())
        }
    }

    struct MockWriter {
        rx: Receiver<String>,
        output: Arc<Mutex<Vec<String>>>,
    }

    impl ChunkWriter for MockWriter {
        type Input = String;

        fn write(self) -> Result<()> {
            for item in self.rx.iter() {
                self.output.lock().unwrap().push(item);
            }
            Ok(())
        }
    }

    impl<'a> ChunkFactory<'a> for MockIO {
        type Reader = MockReader;
        type Parser = MockParser;
        type Writer = MockWriter;

        fn new_reader(&'a self, tx: Sender<String>) -> Result<Self::Reader> {
            Ok(MockReader { tx })
        }

        fn new_parser(
            &'a self,
            rx: Receiver<String>,
            tx: Sender<String>,
        ) -> Result<Self::Parser> {
            Ok(MockParser { rx, tx })
        }

        fn new_writer(&'a self, rx: Receiver<String>) -> Result<Self::Writer> {
            Ok(MockWriter {
                rx,
                output: self.collected.clone(),
            })
        }

        fn channel_reader_parser(
            &'a self,
        ) -> (Sender<String>, Receiver<String>) {
            bounded(4)
        }

        fn channel_parser_writer(
            &'a self,
        ) -> (Sender<String>, Receiver<String>) {
            bounded(4)
        }
    }

    #[test]
    fn test_chunk_io_pipeline() {
        let collected = Arc::new(Mutex::new(Vec::new()));
        let io = MockIO {
            collected: collected.clone(),
        };

        let chunk_io = ChunkIO::new(io);
        chunk_io.run().expect("Pipeline failed");

        let result = collected.lock().unwrap().clone();
        assert_eq!(result, vec!["parsed:line1", "parsed:line2"]);
    }
}
