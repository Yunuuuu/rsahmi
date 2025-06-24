use std::io::{BufWriter, Read, Write};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

use aho_corasick::{AhoCorasick, AhoCorasickKind};
use anyhow::{anyhow, Result};
use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender};
use memchr::{memchr, memrchr};

use super::kractor_match_aho;
use crate::batchsender::BatchSender;

#[allow(clippy::too_many_arguments)]
pub fn reader_kractor_koutput(
    patterns: &[&str],
    file: &str,
    ofile: &str,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    let reader = std::fs::File::open(file)?;
    let matcher = AhoCorasick::builder()
        .kind(Some(AhoCorasickKind::DFA))
        .build(patterns)?;
    let mut writer =
        BufWriter::with_capacity(buffer_size, std::fs::File::create(ofile)?);

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
        let reader = BytesChunkReader::with_capacity(chunk_size, reader);
        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `matcher`
            let has_error = AtomicBool::new(false);
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            rayon::scope(|s| -> Result<()> {
                for chunk_result in reader {
                    let chunk = chunk_result?;
                    if has_error.load(Relaxed) {
                        return Ok(());
                    }
                    s.spawn(|_| {
                        let mut thread_tx = BatchSender::with_capacity(
                            batch_size,
                            parser_tx.clone(),
                        );
                        for line in chunk {
                            if kractor_match_aho(&matcher, &line) {
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

pub struct BytesChunkReader<R>
where
    R: Read,
{
    reader: R,
    chunk_size: usize,
    leftover: BytesMut, // stores bytes after the last \n
}

impl<R> Iterator for BytesChunkReader<R>
where
    R: Read,
{
    type Item = Result<BytesKoutputChunk>;
    fn next(&mut self) -> Option<Self::Item> {
        self.chunk_reader().transpose()
    }
}

impl<R> BytesChunkReader<R>
where
    R: Read,
{
    #[allow(dead_code)]
    pub fn new(reader: R) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Self {
            reader,
            chunk_size: capacity,
            leftover: BytesMut::new(),
        }
    }

    pub fn chunk_reader(&mut self) -> Result<Option<BytesKoutputChunk>> {
        let leftover_len = self.leftover.len();
        let mut buf = BytesMut::with_capacity(self.chunk_size + leftover_len);
        buf.extend_from_slice(&self.leftover);
        self.leftover = BytesMut::new();

        // read files and pass chunks to parser
        unsafe { buf.set_len(leftover_len + self.chunk_size) };
        let nbytes = self.reader.read(&mut buf[leftover_len ..])?;
        unsafe { buf.set_len(leftover_len + nbytes) };
        if nbytes == 0 {
            if buf.is_empty() {
                return Ok(None);
            } else {
                return Ok(Some(BytesKoutputChunk::new(buf.freeze())));
            }
        }

        if let Some(pos) = memrchr(b'\n', &buf) {
            // split at newline
            let chunk = buf.split_to(pos + 1);
            self.leftover = buf; // remainder for next round
            return Ok(Some(BytesKoutputChunk::new(chunk.freeze())));
        } else {
            // no newline found — either final chunk or partial
            // leave all in leftover and try again with larger chunk
            self.chunk_size *= 2;
            self.leftover = buf;
            self.chunk_reader()
        }
    }
}

#[derive(Debug)]
pub struct BytesKoutputChunk {
    pos: usize,
    chunk: Bytes,
}

impl BytesKoutputChunk {
    fn new(chunk: Bytes) -> Self {
        Self { pos: 0, chunk }
    }
}

impl Iterator for BytesKoutputChunk {
    type Item = Bytes;
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.chunk.len() {
            return None;
        }
        let out;
        if let Some(pos) = memchr(b'\n', &self.chunk[self.pos ..]) {
            out = &self.chunk[self.pos ..= self.pos + pos];
            self.pos += pos + 1;
        } else {
            out = &self.chunk[self.pos ..];
            self.pos = self.chunk.len();
        }
        Some(self.chunk.slice_ref(out))
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
        let mut reader = BytesChunkReader::with_capacity(10, cursor);

        let mut all_lines = Vec::new();

        while let Some(chunk_res) = reader.next() {
            let chunk = chunk_res?;
            // Collect all lines from the chunk
            for line in chunk {
                let line_str = std::str::from_utf8(&line).unwrap();
                all_lines.push(line_str.to_string());
            }
        }

        assert_eq!(all_lines, vec![
            "line1", "line2", "line3", "line4", "line5"
        ]);
        Ok(())
    }

    #[test]
    fn test_bytes_chunk_reader_no_final_newline() -> Result<()> {
        let data = b"line1\nline2\nlast_line_without_newline";
        let cursor = Cursor::new(data.as_ref());
        let mut reader = BytesChunkReader::with_capacity(10, cursor);

        let mut all_lines = Vec::new();

        while let Some(chunk_res) = reader.next() {
            let chunk = chunk_res?;
            for line in chunk {
                let line_str = std::str::from_utf8(&line).unwrap();
                all_lines.push(line_str.to_string());
            }
        }

        assert_eq!(all_lines, vec![
            "line1",
            "line2",
            "last_line_without_newline"
        ]);
        Ok(())
    }

    #[test]
    fn test_bytes_chunk_reader_empty() -> Result<()> {
        let data = b"";
        let cursor = Cursor::new(data.as_ref());
        let mut reader = BytesChunkReader::new(cursor);

        match reader.next() {
            None => Ok(()), // expected None means iterator exhausted
            Some(Ok(_)) => {
                Err(anyhow::anyhow!("Expected None but got Some(Ok(_))"))
            }
            Some(Err(e)) => Err(e.into()),
        }
    }
}
