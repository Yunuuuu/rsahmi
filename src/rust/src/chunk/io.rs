use std::io::{Read, Write};

use crossbeam_channel::bounded;

use super::reader::ChunkReader;
use super::ChunkFactory;
use super::{ChunkParser, ChunkProcessor};

pub struct ChunkIO<R, W, F>
where
    R: Read,
    W: Write,
    F: ChunkProcessor,
{
    reader: ChunkReader<R>,
    writer: W,
    factory: ChunkFactory<F>,
    read_queue: usize,
    write_queue: usize,
    threads: usize,
}

impl<R, W, F> ChunkIO<R, W, F>
where
    R: Read,
    W: Write,
    F: ChunkProcessor,
{
    pub fn new(
        reader: ChunkReader<R>,
        writer: W,
        factory: ChunkFactory<F>,
        read_queue: usize,
        write_queue: usize,
        threads: usize,
    ) -> Self {
        Self {
            reader,
            writer,
            factory,
            read_queue,
            write_queue,
            threads,
        }
    }
}

impl<R, W, F> ChunkIO<R, W, F>
where
    R: Read,
    W: Write + Send,
    F: ChunkProcessor,
{
    pub fn process(&mut self) -> std::result::Result<(), String> {
        // two thread is kept for reader and writer
        let threads = if self.threads <= 3 {
            1
        } else {
            self.threads - 2
        };

        let (parser_tx, ref parser_rx) =
            bounded::<Vec<u8>>(threads * self.read_queue);
        let (writer_tx, ref writer_rx) =
            bounded::<Vec<Vec<u8>>>(threads * self.write_queue);

        // Start the processor threads
        std::thread::scope(|scope| {
            // Writer thread: write data to the file
            // accept lines from `writer_rx`
            let writer_handle =
                scope.spawn(|| -> std::result::Result<(), String> {
                    for records in writer_rx.iter() {
                        for record in records {
                            self.writer
                                .write_all(&record)
                                .map_err(|e| format!("Write failed: {e}"))?;
                        }
                    }
                    self.writer
                        .flush()
                        .map_err(|e| format!("Writer flush failed: {e}"))?;
                    Ok(())
                });

            // Worker threads, will send each passed lines to `writer`
            let mut parser_handles = Vec::with_capacity(threads);
            for _ in 0 .. threads {
                let tx = writer_tx.clone();
                let factory = &self.factory;
                let handle =
                    scope.spawn(move || -> std::result::Result<(), String> {
                        let mut buffer = Vec::with_capacity(factory.buffersize);
                        let parser = factory.processor.new_parser();
                        for chunk in parser_rx.iter() {
                            // A list of records is returned
                            parser.parse(chunk, &mut |record| {
                                buffer.push(record);
                                if buffer.len() == buffer.capacity() {
                                    let flushed = std::mem::take(&mut buffer);
                                    tx.send(flushed).map_err(|e| {
                                        format!("Send to writer failed: {}", e)
                                    })?;
                                    buffer =
                                        Vec::with_capacity(factory.buffersize);
                                }
                                Ok(())
                            })?;
                        }
                        if !buffer.is_empty() {
                            tx.send(buffer).map_err(|e| {
                                format!("Flush parsing failed: {}", e)
                            })?;
                        }
                        Ok(())
                    });
                parser_handles.push(handle);
            }

            // read files and pass chunks to parser
            for bytes in &mut self.reader {
                match bytes {
                    Ok(chunk) => {
                        // send chunk to parser
                        parser_tx.send(chunk).map_err(|e| {
                            format!("Send to parser failed: {e}")
                        })?;
                    }
                    Err(e) => return Err(e),
                }
            }
            // close parser channel to stop parsers
            drop(parser_tx);

            // we ensure no data to send
            for handler in parser_handles {
                let _ = handler.join().map_err(|_| "parser thread panicked")?;
            }
            drop(writer_tx); // close writer channel to stop writer

            // we ensure all lines have been writen
            let _ =
                writer_handle.join().map_err(|_| "writer thread panicked")?;
            Ok(())
        })
    }
}
