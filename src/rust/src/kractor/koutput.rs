use std::fs::File;
use std::io::{BufWriter, Write};

use aho_corasick::AhoCorasick;
use anyhow::{anyhow, Ok, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use memchr::{memchr, memrchr};

use super::buffer::Buffer;
use super::io::*;
use crate::kractor::batchsender::BatchSender;

#[allow(clippy::too_many_arguments)]
pub fn kractor_koutput(
    file: &str,
    ofile: &str,
    matcher: &AhoCorasick,
    read_buffer: usize,
    write_buffer: usize,
    batch_size: usize,
    read_queue: usize,
    write_queue: usize,
    threads: usize,
) -> Result<()> {
    ChunkIO::new(KoutputFactory::new(
        file,
        ofile,
        matcher,
        read_buffer,
        write_buffer,
        batch_size,
        read_queue,
        write_queue,
        threads,
    ))
    .run()
}

struct KoutputFactory<'a> {
    file: &'a str,
    ofile: &'a str,
    matcher: &'a AhoCorasick,
    read_buffer: usize,
    write_buffer: usize,
    batch_size: usize,
    read_queue: usize,
    write_queue: usize,
    threads: usize,
}

impl<'a> KoutputFactory<'a> {
    fn new(
        file: &'a str,
        ofile: &'a str,
        matcher: &'a AhoCorasick,
        read_buffer: usize,
        write_buffer: usize,
        batch_size: usize,
        read_queue: usize,
        write_queue: usize,
        threads: usize,
    ) -> Self {
        // two thread is kept for reader and writer
        let threads = if threads <= 3 { 1 } else { threads - 2 };
        Self {
            file,
            ofile,
            matcher,
            read_buffer,
            write_buffer,
            batch_size,
            read_queue,
            write_queue,
            threads,
        }
    }
}

impl<'a> ChunkFactory<'a> for KoutputFactory<'a> {
    type Reader = KoutputReader;
    type Parser = KoutputParser<'a>;
    type Writer = KoutputWriter;
    fn new_reader(
        &'a self,
        tx: Sender<<Self::Reader as ChunkReader>::Output>,
    ) -> Result<Self::Reader> {
        let reader = std::fs::File::open(self.file)?;
        let buffer = Buffer::with_capacity(self.read_buffer);
        Ok(KoutputReader {
            reader,
            buf: buffer,
            capacity: self.read_buffer,
            tx: tx,
        })
    }

    fn new_parser(
        &'a self,
        rx: Receiver<<Self::Parser as ChunkParser>::Input>,
        tx: Sender<<Self::Parser as ChunkParser>::Output>,
    ) -> Result<Self::Parser> {
        let tx = BatchSender::with_capacity(self.batch_size, tx);
        Ok(KoutputParser {
            rx: rx,
            tx: tx,
            matcher: &self.matcher,
            threads: self.threads,
        })
    }

    fn new_writer(
        &'a self,
        rx: Receiver<<Self::Writer as ChunkWriter>::Input>,
    ) -> Result<Self::Writer> {
        let writer = std::fs::File::create(self.ofile)?;
        Ok(KoutputWriter {
            rx,
            writer: BufWriter::with_capacity(self.write_buffer, writer),
        })
    }

    fn channel_reader_parser(
        &'a self,
    ) -> (
        Sender<<Self::Reader as ChunkReader>::Output>,
        Receiver<<Self::Reader as ChunkReader>::Output>,
    ) {
        bounded(self.read_queue * self.threads)
    }

    fn channel_parser_writer(
        &'a self,
    ) -> (
        Sender<<Self::Parser as ChunkParser>::Output>,
        Receiver<<Self::Parser as ChunkParser>::Output>,
    ) {
        bounded(self.write_queue * self.threads)
    }
}

pub struct KoutputReader {
    reader: File,
    buf: Buffer,
    capacity: usize,
    tx: Sender<Vec<u8>>,
}

impl KoutputReader {
    fn take_chunk(&mut self) -> Result<Option<Vec<u8>>> {
        self.buf.compact();
        loop {
            match self.buf.fill_buf(&mut self.reader)? {
                Some(_) => {
                    let buf = self.buf.buffer();
                    if let Some(pos) = memrchr(b'\n', buf) {
                        let chunk = buf[..= pos].to_vec(); // including the final '\n'

                        self.consume(chunk.len());
                        return Ok(Some(chunk));
                    } else {
                        self.buf.extend(self.capacity);
                    }
                }
                None => {
                    // end of the file
                    if self.buf.filled() == 0 {
                        // empty buffer
                        return Ok(None);
                    } else {
                        // return all remaining buffer
                        let chunk = self.buf.buffer().to_vec();
                        self.consume(chunk.len());
                        return Ok(Some(chunk));
                    }
                }
            }
        }
    }

    fn send(&mut self, msg: Vec<u8>) -> Result<()> {
        self.tx.send(msg)?;
        Ok(())
    }

    fn consume(&mut self, amt: usize) {
        self.buf.consume(amt);
    }
}

impl ChunkReader for KoutputReader {
    type Output = Vec<u8>;
    fn read(mut self) -> Result<()> {
        // read files and pass chunks to parser
        while let Some(chunk) = self.take_chunk()? {
            self.send(chunk)?;
        }
        Ok(())
    }
}

pub struct KoutputParser<'a> {
    rx: Receiver<Vec<u8>>,
    tx: BatchSender<Vec<u8>>,
    matcher: &'a AhoCorasick,
    threads: usize,
}

impl<'a> ChunkParser for KoutputParser<'a> {
    type Input = Vec<u8>;
    type Output = Vec<Vec<u8>>;
    fn parse(mut self) -> Result<()> {
        if self.threads <= 1 {
            for chunk in &self.rx {
                for line in KoutputMatcherLines::new(self.matcher, chunk) {
                    self.tx.send(line)?;
                }
            }
            self.tx.flush()?;
            Ok(())
        } else {
            std::thread::scope(|scope| {
                let mut handles = Vec::with_capacity(self.threads);
                for _ in 0 .. self.threads {
                    let handle = scope.spawn(|| -> Result<()> {
                        let mut tx = self.tx.clone();
                        for chunk in &self.rx {
                            for line in
                                KoutputMatcherLines::new(self.matcher, chunk)
                            {
                                tx.send(line)?;
                            }
                        }
                        tx.flush()?;
                        Ok(())
                    });
                    handles.push(handle);
                }
                for handle in handles {
                    handle.join().map_err(|e| {
                        anyhow!("Parser scoped thread panicked: {:?}", e)
                    })??;
                }
                Ok(())
            })
        }
    }
}

struct KoutputMatcherLines<'a> {
    matcher: &'a AhoCorasick,
    chunk: Vec<u8>,
    position: usize,
}

impl<'a> KoutputMatcherLines<'a> {
    fn new(matcher: &'a AhoCorasick, chunk: Vec<u8>) -> Self {
        Self {
            matcher,
            chunk,
            position: 0,
        }
    }

    fn match_line(&self, line: &[u8]) -> bool {
        // Efficient 3rd column parsing
        let mut field_start = 0usize;
        let mut field_count = 0usize;
        while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
            if field_count == 2 {
                // we don't include the last `\t`
                let taxid = &line[field_start .. (field_start + tab_pos)];
                return self.matcher.find(taxid).is_some();
            }
            field_start += tab_pos + 1;
            field_count += 1;
        }
        false
    }
}

impl<'a> Iterator for KoutputMatcherLines<'a> {
    type Item = Vec<u8>; // underscore lifetime means "tied to &mut self"

    fn next(&mut self) -> Option<Self::Item> {
        while self.position < self.chunk.len() {
            let start = self.position;
            match memchr(b'\n', &self.chunk[start ..]) {
                Some(pos) => {
                    self.position += pos + 1;
                }
                None => {
                    self.position = self.chunk.len();
                }
            };
            let line = &self.chunk[start .. self.position];
            if self.match_line(line) {
                return Some(line.to_vec());
            }
        }
        None
    }
}

pub struct KoutputWriter {
    rx: Receiver<Vec<Vec<u8>>>,
    writer: BufWriter<File>,
}

impl ChunkWriter for KoutputWriter {
    type Input = Vec<Vec<u8>>;
    fn write(self) -> Result<()> {
        let mut writer = self.writer;
        for batch in self.rx.iter() {
            for record in batch {
                writer.write_all(&record)?;
            }
        }
        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::fs::{remove_file, File};
    use std::io::{Read, Write};

    use aho_corasick::AhoCorasick;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_kractor_koutput_pipeline() -> Result<()> {
        // Step 1: Create test input file
        let mut input_file = NamedTempFile::new()?;
        writeln!(input_file, "a\tb\t12345\tmore")?;
        writeln!(input_file, "a\tb\t67890\tmore")?;
        writeln!(input_file, "a\tb\t11111\tmore")?;
        writeln!(input_file, "a\tb\t99999\tmore")?;
        let input_path = input_file.into_temp_path();

        // Step 2: Create output file path
        let output_file = NamedTempFile::new()?;
        let output_path = output_file.into_temp_path();

        // Step 3: Aho-Corasick matcher for 12345 and 11111
        let matcher = AhoCorasick::new(["12345", "11111"])?;

        // Step 4: Run the Koutput pipeline
        kractor_koutput(
            input_path.to_str().unwrap(),
            output_path.to_str().unwrap(),
            &matcher,
            64, // read buffer
            64, // write buffer
            2,  // batch size
            4,  // read queue
            4,  // write queue
            4,  // threads
        )?;

        // Step 5: Read output and verify
        let mut output = String::new();
        File::open(&output_path)?.read_to_string(&mut output)?;

        let lines: Vec<&str> = output.trim().lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines.iter().any(|line| line.contains("12345")));
        assert!(lines.iter().any(|line| line.contains("11111")));
        assert!(!lines.iter().any(|line| line.contains("67890")));
        assert!(!lines.iter().any(|line| line.contains("99999")));

        // Cleanup
        let _ = remove_file(input_path);
        let _ = remove_file(output_path);

        Ok(())
    }
}
