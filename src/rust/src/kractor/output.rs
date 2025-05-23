use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use aho_corasick::AhoCorasick;
use crossbeam_channel::bounded;
use extendr_api::prelude::*;
use memchr::memchr;

#[allow(clippy::too_many_arguments)]
pub fn write_matching_output<P>(
    koutput: P,
    matcher: &AhoCorasick,
    ofile: P,
    io_buffer: usize,
    buffersize: usize,
    batchsize: usize,
    queue_capacity: usize,
    threads: usize,
) -> std::result::Result<(), String>
where
    P: AsRef<Path> + Display,
{
    rprintln!("Extracting matching kraken2 output from: {}", koutput);
    // one thread is kept for writer
    let threads = if threads <= 2 { 1 } else { threads - 1 };
    let input = File::open(koutput).map_err(|e| e.to_string())?;
    let mut output = BufWriter::with_capacity(
        io_buffer,
        File::create(ofile).map_err(|e| e.to_string())?,
    );

    // Start the processor threads
    let (writer_tx, ref writer_rx) =
        bounded::<Vec<Vec<u8>>>(threads * queue_capacity);
    let (work_tx, ref work_rx) = bounded::<Vec<u8>>(threads * queue_capacity);

    std::thread::scope(|scope| {
        // let patterns = Arc::new(patterns);

        // Writer thread: write data to the file
        // accept lines from `writer_rx`
        let writer_handle =
            scope.spawn(move || -> std::result::Result<(), String> {
                for batch in writer_rx.iter() {
                    for line in batch {
                        output
                            .write_all(&line)
                            .map_err(|e| format!("Write failed: {e}"))?;
                    }
                }
                output.flush().map_err(|e| format!("Flush failed: {e}"))?;
                Ok(())
            });

        // Worker threads, will send each passed lines to `writer`
        let mut worker_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let tx = writer_tx.clone();
            let handle =
                scope.spawn(move || -> std::result::Result<(), String> {
                    let mut processor = KOutputChunkParser::new(batchsize, tx);
                    for chunk in work_rx.iter() {
                        // we don't include the last `\n`
                        processor.import_chunk(&chunk, matcher)?;
                    }
                    processor.close()
                });
            worker_handles.push(handle);
        }

        // read files and pass lines
        let mut reader = BytesReader::new(buffersize, input, work_tx.clone());
        loop {
            let read_bytes = reader.read_buffer()?;
            if read_bytes == 0 {
                reader.close()?;
                break;
            }
            reader.send_buffer()?;
        }

        // close work channel to stop workers
        drop(work_tx);

        // we ensure no lines to send
        for handler in worker_handles {
            let _ = handler.join().map_err(|_| "worker thread panicked")?;
        }
        drop(writer_tx); // close writer channel to stop writer

        // we ensure all lines have been writen
        let _ = writer_handle.join().map_err(|_| "writer thread panicked")?;
        Ok(())
    })
}

struct KOutputChunkParser {
    // size of the buffer to read
    buffersize: usize,
    // buffer to store lines
    buffer: Vec<Vec<u8>>,
    // used to send lines to writer
    // we use a channel to avoid blocking the worker threads
    sender: crossbeam_channel::Sender<Vec<Vec<u8>>>,
}

impl KOutputChunkParser {
    fn new(
        buffersize: usize,
        sender: crossbeam_channel::Sender<Vec<Vec<u8>>>,
    ) -> Self {
        Self {
            buffersize,
            buffer: Vec::with_capacity(buffersize),
            sender,
        }
    }

    fn close(self) -> std::result::Result<(), String> {
        if !self.buffer.is_empty() {
            self.sender
                .send(self.buffer)
                .map_err(|e| format!("Flush failed: {e}"))?;
        }
        Ok(())
    }

    fn send(&mut self) -> std::result::Result<(), String> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let mut out = Vec::with_capacity(self.buffer.capacity());
        std::mem::swap(&mut out, &mut self.buffer);
        self.sender
            .send(out)
            .map_err(|e| format!("Send to writer failed: {e}"))
    }

    fn push(&mut self, line: Vec<u8>) -> std::result::Result<(), String> {
        if self.buffer.len() >= self.buffersize {
            self.send()?;
        }
        self.buffer.push(line);
        Ok(())
    }

    fn import_chunk(
        &mut self,
        chunk: &[u8],
        matcher: &AhoCorasick,
    ) -> std::result::Result<(), String> {
        let mut start = 0;
        while let Some(line_pos) = memchr(b'\n', &chunk[start ..]) {
            // we include the last `\n` for writing
            let line = &chunk[start ..= start + line_pos];
            self.import_line(line, matcher)?;
            start += line_pos + 1;
        }
        Ok(())
    }

    fn import_line(
        &mut self,
        line: &[u8],
        matcher: &AhoCorasick,
    ) -> std::result::Result<(), String> {
        // Efficient 3rd column parsing
        let mut field_start = 0usize;
        let mut field_count = 0usize;
        while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
            if field_count == 2 {
                // we don't include the last `\t`
                let taxid = &line[field_start .. field_start + tab_pos];
                if matcher.find(taxid).is_some() {
                    self.push(line.to_vec())?;
                }
                break;
            }
            field_start = tab_pos + 1;
            field_count += 1;
        }
        Ok(())
    }
}

struct BytesReader<R>
where
    R: Read,
{
    reader: R,
    // size of the buffer to read
    buffersize: usize,
    // buffer to store lines
    buffer: Vec<u8>,
    // used to store the byte from where to read
    offset: usize,
    // used to send bytes
    channel: crossbeam_channel::Sender<Vec<u8>>,
}

impl<R> BytesReader<R>
where
    R: Read,
{
    fn new(
        buffersize: usize,
        reader: R,
        channel: crossbeam_channel::Sender<Vec<u8>>,
    ) -> Self {
        Self {
            buffersize,
            buffer: vec![0; buffersize],
            offset: 0,
            reader,
            channel,
        }
    }

    fn read_buffer(&mut self) -> std::result::Result<usize, String> {
        let read_bytes = self
            .reader
            .read(&mut self.buffer[self.offset ..])
            .map_err(|e| format!("Read failed: {e}"))?;
        self.offset += read_bytes;
        Ok(read_bytes)
    }

    fn send_buffer(&mut self) -> std::result::Result<(), String> {
        match self.buffer[.. self.offset]
            .iter()
            .rposition(|b| *b == b'\n')
        {
            Some(pos) => {
                let chunk = self.buffer.drain(..= pos).collect::<Vec<u8>>();
                self.channel
                    .send(chunk)
                    .map_err(|e| format!("Send to workers failed: {e}"))?;
                // drain won't reduce the capacity of the buffer but will remove the data
                self.fill_buffer();
                // we need to move the offset to the end of the remaing buffer
                self.offset -= pos + 1;
            }
            None => {
                if self.offset == self.buffer.capacity() {
                    self.extend_buffer()?;
                }
            }
        };
        Ok(())
    }

    fn close(&mut self) -> std::result::Result<(), String> {
        // ensure all buffer has been send
        if self.offset > 0 {
            let mut chunk =
                self.buffer.drain(..= self.offset).collect::<Vec<u8>>();
            // always ensure the last line has `\n`
            if !chunk.ends_with(b"\n") {
                chunk.push(b'\n');
            }
            self.channel
                .send(chunk)
                .map_err(|e| format!("Send to workers failed: {e}"))?;
        }
        Ok(())
    }

    // everytime we increase the buffer, we need to fill it with 0
    fn fill_buffer(&mut self) {
        self.buffer.resize(self.buffer.capacity(), 0);
    }

    fn extend_buffer(&mut self) -> std::result::Result<(), String> {
        self.buffer
            .try_reserve(self.buffersize)
            .map_err(|e| format!("Increase buffer failed: {e}"))?;
        self.fill_buffer();
        Ok(())
    }
}
