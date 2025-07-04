use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

use aho_corasick::AhoCorasick;
use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use memchr::memrchr;

use super::kractor_match_aho;
use crate::batchsender::BatchSender;
use crate::reader::bytes::{BytesLineReader, BytesProgressBarReader};

#[allow(clippy::too_many_arguments)]
pub fn reader_kractor_koutput<R: Read + Send>(
    include_aho: AhoCorasick,
    exclude_aho: Option<AhoCorasick>,
    reader: BytesProgressBarReader<R>,
    ofile: &str,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    let mut writer = BufWriter::with_capacity(buffer_size, File::create(ofile)?);

    std::thread::scope(|scope| {
        // Create a channel between the parser and writer threads
        // The channel transmits batches
        let (parser_tx, writer_rx): (Sender<Vec<Bytes>>, Receiver<Vec<Bytes>>) =
            crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        let writer_handle = scope.spawn(move || -> Result<()> {
            for chunk in writer_rx {
                for line in chunk {
                    writer.write_all(&line)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let mut reader = KoutputBytesChunkReader::with_capacity(chunk_size, reader);
        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, `include_aho` and `exclude_aho`
            let has_error = AtomicBool::new(false);
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            rayon::scope(|s| -> Result<()> {
                while let Some(chunk) = reader.chunk_reader()? {
                    if has_error.load(Relaxed) {
                        return Ok(());
                    }
                    s.spawn(|_| {
                        let mut chunk = chunk; // move the chunk
                        let mut thread_tx =
                            BatchSender::with_capacity(batch_size, parser_tx.clone());
                        while let Some(line) = chunk.read_line_inclusive() {
                            if kractor_match_aho(&include_aho, &exclude_aho, &line) {
                                match thread_tx.send(line) {
                                    Ok(_) => continue,
                                    Err(e) => {
                                        has_error.store(true, Relaxed);
                                        // we only capture the first error
                                        let _ = err_tx.try_send(e);
                                        return ();
                                    }
                                };
                            }
                        }
                        match thread_tx.flush() {
                            Ok(_) => return (),
                            Err(e) => {
                                has_error.store(true, Relaxed);
                                // we only capture the first error
                                let _ = err_tx.try_send(e);
                                return ();
                            }
                        }
                    });
                }
                Ok(())
            })?;
            drop(err_tx);
            drop(parser_tx);
            if has_error.load(Relaxed) {
                let err = err_rx.recv()?;
                Err(anyhow!("Failed to send to Writer thread: {}", err))
            } else {
                Ok(())
            }
        });

        // ─── Join Threads and Propagate Errors ────────────────
        writer_handle
            .join()
            .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))??;
        parser_handle
            .join()
            .map_err(|e| anyhow!("Parser thread panicked: {:?}", e))??;
        Ok(())
    })
}

pub struct KoutputBytesChunkReader<R>
where
    R: Read,
{
    chunk_size: usize,
    reader: BytesProgressBarReader<R>,
}

impl<R: Read> KoutputBytesChunkReader<R> {
    #[allow(dead_code)]
    pub fn new(reader: BytesProgressBarReader<R>) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub fn with_capacity(capacity: usize, reader: BytesProgressBarReader<R>) -> Self {
        Self {
            chunk_size: capacity,
            reader,
        }
    }

    pub fn chunk_reader(&mut self) -> Result<Option<BytesLineReader>> {
        let buf = match self.reader.read_bytes(self.chunk_size)? {
            Some(data) => data,
            None => return Ok(None),
        };

        let chunk;
        if buf.eof() {
            chunk = buf.into_bytes();
        } else {
            let mut bytes = buf.into_bytes();
            if let Some(pos) = memrchr(b'\n', &bytes) {
                chunk = bytes.split_to(pos + 1);
                self.reader.take_leftover(bytes);
            } else {
                // no newline found — either final chunk or partial
                // leave all in leftover and try again with larger chunk
                self.reader.take_leftover(bytes);
                self.chunk_size *= 2;
                return self.chunk_reader();
            }
        }
        return Ok(Some(BytesLineReader::new(chunk.freeze())));
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use anyhow::Result;

    use super::*;

    #[test]
    fn test_bytes_chunk_reader_basic() -> Result<()> {
        let data = b"line1\nline2\nline3\nline4\nline5\n";
        let cursor = Cursor::new(data.as_ref());
        let mut reader =
            KoutputBytesChunkReader::with_capacity(10, BytesProgressBarReader::new(cursor));

        let mut all_lines = Vec::new();

        while let Some(mut chunk) = reader.chunk_reader()? {
            // Collect all lines from the chunk
            while let Some(line) = chunk.read_line_inclusive() {
                // println!("Reading line {}", String::from_utf8_lossy(&line));
                let line_str = std::str::from_utf8(&line).unwrap();
                all_lines.push(line_str.to_string());
            }
        }

        assert_eq!(all_lines, vec![
            "line1\n", "line2\n", "line3\n", "line4\n", "line5\n"
        ]);
        Ok(())
    }

    #[test]
    fn test_bytes_chunk_reader_no_final_newline() -> Result<()> {
        let data = b"line1\nline2\nlast_line_without_newline";
        let cursor = Cursor::new(data.as_ref());
        let mut reader =
            KoutputBytesChunkReader::with_capacity(10, BytesProgressBarReader::new(cursor));

        let mut all_lines = Vec::new();

        while let Some(mut chunk) = reader.chunk_reader()? {
            // println!("Reading bytes: {}", String::from_utf8_lossy(chunk.borrow_bytes()));
            while let Some(line) = chunk.read_line_inclusive() {
                // println!("Reading line {}", String::from_utf8_lossy(&line));
                let line_str = std::str::from_utf8(&line).unwrap();
                all_lines.push(line_str.to_string());
            }
        }

        assert_eq!(all_lines, vec![
            "line1\n",
            "line2\n",
            "last_line_without_newline"
        ]);
        Ok(())
    }

    #[test]
    fn test_bytes_chunk_reader_empty() -> Result<()> {
        let data = b"";
        let cursor = Cursor::new(data.as_ref());
        let mut reader = KoutputBytesChunkReader::new(BytesProgressBarReader::new(cursor));

        match reader.chunk_reader()? {
            None => Ok(()), // expected None means iterator exhausted
            Some(_) => Err(anyhow::anyhow!("Expected None but got Some(Ok(_))")),
        }
    }
}
