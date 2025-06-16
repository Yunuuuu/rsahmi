use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::sync::atomic::{AtomicBool, Ordering};

use aho_corasick::AhoCorasick;
use anyhow::{anyhow, Result};
use bytes::{Bytes, BytesMut};
use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use memchr::memchr;
use memmap2::{Advice, Mmap};

use super::io::*;
use crate::kractor::batchsender::BatchSender;

fn new_channel<T>(queue: Option<usize>) -> (Sender<T>, Receiver<T>) {
    if let Some(queue) = queue {
        bounded(queue)
    } else {
        unbounded()
    }
}

#[allow(clippy::too_many_arguments)]
pub fn kractor_koutput(
    file: &str,
    ofile: &str,
    matcher: &AhoCorasick,
    read_buffer: usize,
    write_buffer: usize,
    batch_size: usize,
    read_queue: Option<usize>,
    write_queue: Option<usize>,
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
    ))
    .run()
}

pub fn mmap_kractor_koutput(
    file: &str,
    ofile: &str,
    matcher: &AhoCorasick,
    batch_size: usize,
    write_buffer: usize,
    read_queue: Option<usize>,
    write_queue: Option<usize>,
) -> Result<()> {
    let file = File::open(file)?;
    // https://github.com/rayon-rs/rayon/discussions/1164
    let map = unsafe { Mmap::map(&file) }?;
    map.advise(Advice::Sequential)?;
    let mut writer =
        BufWriter::with_capacity(write_buffer, std::fs::File::create(ofile)?);
    std::thread::scope(|scope| {
        let (reader_tx, parser_rx) = new_channel(read_queue);
        let (parser_tx, writer_rx) = new_channel(write_queue);
        let writer_handle = scope.spawn(|| -> Result<()> {
            for chunk in writer_rx {
                for line in chunk {
                    writer.write_all(line)?;
                }
            }
            Ok(())
        });
        let parser_handle = scope.spawn(|| -> Result<()> {
            let parser_err = AtomicBool::new(false); // Shared error flag
            rayon::scope(|s| {
                let parser_batch_tx =
                    BatchSender::with_capacity(batch_size, parser_tx);
                for chunk in parser_rx {
                    let mut tx = parser_batch_tx.clone();
                    let parser_err_ref = &parser_err;
                    s.spawn(move |_| {
                        for line in chunk {
                            if kractor_match_line(matcher, line) {
                                match tx.send(line) {
                                    Err(_) => {
                                        if parser_err_ref
                                            .load(Ordering::Relaxed)
                                        {
                                            return ();
                                        } else {
                                            parser_err_ref
                                                .store(true, Ordering::Relaxed)
                                        }
                                    }
                                    Ok(_) => continue,
                                };
                            }
                        }
                        if let Err(_) = tx.flush() {
                            if !parser_err_ref.load(Ordering::Relaxed) {
                                parser_err_ref.store(true, Ordering::Relaxed)
                            }
                        }
                    });
                }
            });
            if parser_err.load(Ordering::Relaxed) {
                Err(anyhow!("Failed to send to Writer thread"))
            } else {
                Ok(())
            }
        });
        let reader_handle = scope.spawn(|| -> Result<()> {
            let mut start = 0usize;
            let mut reader_batch_tx =
                BatchSender::with_capacity(batch_size, reader_tx);
            while let Some(pos) = memchr(b'\n', &map[start ..]) {
                let line = &map[start ..= start + pos];
                reader_batch_tx.send(line).map_err(|e| {
                    anyhow!("Failed to send to Parser thread: {}", e)
                })?;
                start += pos + 1
            }
            let line = &map[start ..];
            reader_batch_tx.send(line).map_err(|e| {
                anyhow!("Failed to send to Parser thread: {}", e)
            })?;
            reader_batch_tx.flush().map_err(|e| {
                anyhow!("Failed to send to Parser thread: {}", e)
            })?;
            Ok(())
        });
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

#[cfg(feature = "bench")]
pub fn bench_kractor_koutput_reader(
    file: &str,
    buffer_size: usize,
    n_queue: Option<usize>,
) -> Result<()> {
    use std::time::Instant;

    use anyhow::anyhow;

    let file = std::fs::File::open(file)?;
    let (tx, rx) = new_channel(n_queue);
    let mut kreader = KoutputReader {
        reader: file,
        capacity: buffer_size,
        tx: tx,
    };
    let handle = std::thread::spawn(move || {
        let mut count = 0;
        for chunk in rx {
            count += chunk.len()
        }
        count
    });
    let start = Instant::now();
    kreader.read()?;
    drop(kreader);
    let nbytes = handle
        .join()
        .map_err(|e| anyhow!("Thread panicked: {:?}", e))?;
    let elapsed = start.elapsed();
    println!(
        "Total bytes: {} in {:.2?} → {:.2} MB/s | {:.2} GB/s",
        nbytes,
        elapsed,
        nbytes as f64 / 1024.0 / 1024.0 / elapsed.as_secs_f64(),
        nbytes as f64 / 1024.0 / 1024.0 / 1024.0 / elapsed.as_secs_f64()
    );
    Ok(())
}

struct KoutputFactory<'a> {
    file: &'a str,
    ofile: &'a str,
    matcher: &'a AhoCorasick,
    read_buffer: usize,
    write_buffer: usize,
    batch_size: usize,
    read_queue: Option<usize>,
    write_queue: Option<usize>,
}

impl<'a> KoutputFactory<'a> {
    fn new(
        file: &'a str,
        ofile: &'a str,
        matcher: &'a AhoCorasick,
        read_buffer: usize,
        write_buffer: usize,
        batch_size: usize,
        read_queue: Option<usize>,
        write_queue: Option<usize>,
    ) -> Self {
        Self {
            file,
            ofile,
            matcher,
            read_buffer,
            write_buffer,
            batch_size,
            read_queue,
            write_queue,
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
        Ok(KoutputReader {
            reader,
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
        new_channel(self.read_queue)
    }

    fn channel_parser_writer(
        &'a self,
    ) -> (
        Sender<<Self::Parser as ChunkParser>::Output>,
        Receiver<<Self::Parser as ChunkParser>::Output>,
    ) {
        new_channel(self.write_queue)
    }
}

pub struct KoutputReader {
    reader: File,
    capacity: usize,
    tx: Sender<Bytes>,
}

impl ChunkReader for KoutputReader {
    type Output = Bytes;

    #[inline]
    fn read(&mut self) -> Result<()> {
        let mut buf;
        loop {
            buf = BytesMut::with_capacity(self.capacity);
            // read files and pass chunks to parser
            unsafe { buf.set_len(self.capacity) };
            let nbytes = self.reader.read(&mut buf)?;
            if nbytes == 0 {
                break;
            }
            if nbytes < buf.len() {
                // SAFETY: Shrinking the buffer cannot expose uninitialized bytes.
                unsafe { buf.set_len(nbytes) };
            }
            self.tx.send(buf.freeze())?;
        }
        Ok(())
    }
}

pub struct KoutputParser<'a> {
    rx: Receiver<Bytes>,
    tx: BatchSender<Vec<u8>>,
    matcher: &'a AhoCorasick,
}

impl<'a> KoutputParser<'a> {
    fn send_line_slice(&mut self, line: &[u8]) -> Result<()> {
        if self.match_line(line) {
            self.send_slice(line)?;
        }
        Ok(())
    }
    fn send_line_vec(&mut self, line: Vec<u8>) -> Result<()> {
        if self.match_line(&line) {
            self.send_vec(line)?;
        }
        Ok(())
    }
    fn flush(&mut self) -> Result<()> {
        self.tx.flush()?;
        Ok(())
    }
    fn send_slice(&mut self, line: &[u8]) -> Result<()> {
        self.send_vec(line.to_vec())
    }
    fn send_vec(&mut self, line: Vec<u8>) -> Result<()> {
        self.tx.send(line)?;
        Ok(())
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

impl<'a> ChunkParser for KoutputParser<'a> {
    type Input = Bytes;
    type Output = Vec<Vec<u8>>;
    fn parse(&mut self) -> Result<()> {
        let mut leftover: Option<Bytes> = None;
        let mut start;
        while let std::result::Result::Ok(chunk) = self.rx.recv() {
            if let Some(prev) = leftover {
                // Combine previous partial line with new chunk — 1 allocation, unavoidable
                if let Some(pos) = memchr(b'\n', &chunk) {
                    let line: Vec<u8> = prev
                        .iter()
                        .chain(chunk.slice(..= pos).iter())
                        .copied()
                        .collect();
                    self.send_line_vec(line)?;
                    start = pos + 1;
                } else {
                    let line: Vec<u8> =
                        prev.iter().chain(chunk.iter()).copied().collect();
                    leftover = Some(Bytes::from(line));
                    continue;
                }
            } else {
                start = 0;
            }
            while let Some(pos) = memchr(b'\n', &chunk[start ..]) {
                let end = start + pos + 1;
                self.send_line_slice(&chunk[start .. end])?;
                start = end;
            }
            if start < chunk.len() {
                leftover = Some(chunk.slice(start ..)); // save remainder (no copy)
            } else {
                leftover = None
            }
        }
        if let Some(line) = leftover {
            self.send_line_vec(Vec::from(line))?;
        }
        self.flush()?;
        Ok(())
    }
}

pub struct KoutputWriter {
    rx: Receiver<Vec<Vec<u8>>>,
    writer: BufWriter<File>,
}

impl ChunkWriter for KoutputWriter {
    type Input = Vec<Vec<u8>>;
    fn write(&mut self) -> Result<()> {
        for batch in self.rx.iter() {
            for record in batch {
                self.writer.write_all(&record)?;
            }
        }
        self.writer.flush()?;
        Ok(())
    }
}

fn kractor_match_line(matcher: &AhoCorasick, line: &[u8]) -> bool {
    // Efficient 3rd column parsing
    let mut field_start = 0usize;
    let mut field_count = 0usize;
    while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
        if field_count == 2 {
            // we don't include the last `\t`
            let taxid = &line[field_start .. (field_start + tab_pos)];
            return matcher.find(taxid).is_some();
        }
        field_start += tab_pos + 1;
        field_count += 1;
    }
    false
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
            64,      // read buffer
            64,      // write buffer
            2,       // batch size
            Some(4), // read queue
            Some(4), // write queue
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

    #[test]
    #[cfg(feature = "bench")]
    fn bench_bench_kractor_koutput_reader() -> Result<()> {
        bench_kractor_koutput_reader(
            "../../bench/data/CNP000460_P01N_output.txt",
            1 * 1024 * 1024,
            None,
        )
    }
}
