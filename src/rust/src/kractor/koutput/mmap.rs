use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

use aho_corasick::{AhoCorasick, AhoCorasickKind};
use anyhow::{anyhow, Result};
use memchr::{memchr, memrchr};
use memmap2::{Advice, Mmap};

use super::kractor_match_aho;
use crate::batchsender::BatchSender;

pub fn mmap_kractor_koutput(
    patterns: &[&str],
    file: &str,
    ofile: &str,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    let file = File::open(file)?;
    let matcher = AhoCorasick::builder()
        .kind(Some(AhoCorasickKind::DFA))
        .build(patterns)?;

    let mut writer =
        BufWriter::with_capacity(buffer_size, std::fs::File::create(ofile)?);

    // let taxid_sets: HashSet<&[u8]> =
    //     patterns.into_iter().map(|s| s.as_bytes()).collect();
    // https://github.com/rayon-rs/rayon/discussions/1164
    let map = unsafe { Mmap::map(&file) }?;
    map.advise(Advice::Sequential)?;

    std::thread::scope(|scope| {
        // Create a channel between the parser and writer threads
        // The channel transmits batches
        let (parser_tx, writer_rx) = crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        let writer_handle = scope.spawn(move || -> Result<()> {
            for chunk in writer_rx {
                for line in chunk {
                    writer.write_all(line)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let reader = SliceChunkReader::with_capacity(chunk_size, &map);
        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `matcher`
            let has_error = AtomicBool::new(false);
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            rayon::scope(|s| -> Result<()> {
                for chunk in reader {
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
                            Err(_) => {
                                has_error.store(true, Relaxed);
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

pub struct SliceChunkReader<'a> {
    pos: usize,
    slice: &'a [u8],
    chunk_size: usize,
}

impl<'a> Iterator for SliceChunkReader<'a> {
    type Item = SliceKoutputChunk<'a>;
    fn next(&mut self) -> Option<Self::Item> {
        self.chunk_reader()
    }
}

impl<'a> SliceChunkReader<'a> {
    #[allow(dead_code)]
    pub fn new(slice: &'a [u8]) -> Self {
        Self::with_capacity(8 * 1024, slice)
    }

    pub fn with_capacity(capacity: usize, slice: &'a [u8]) -> Self {
        Self {
            pos: 0,
            slice,
            chunk_size: capacity,
        }
    }

    pub fn chunk_reader(&mut self) -> Option<SliceKoutputChunk<'a>> {
        // If we've reached or passed the end of the input buffer, stop reading
        if self.pos >= self.slice.len() {
            return None;
        }
        let end = self.pos + self.chunk_size;
        let mut chunk;
        if end >= self.slice.len() {
            chunk = &self.slice[self.pos ..];
            self.pos = self.slice.len();
        } else {
            chunk = &self.slice[self.pos .. end];
            if let Some(pos) = memrchr(b'\n', chunk) {
                chunk = &chunk[..= pos];
                self.pos += pos + 1;
            } else {
                self.chunk_size *= 2;
                return self.chunk_reader();
            }
        }
        Some(SliceKoutputChunk::new(chunk))
    }
}

#[derive(Debug)]
pub struct SliceKoutputChunk<'a> {
    pos: usize,
    chunk: &'a [u8],
}

impl<'a> SliceKoutputChunk<'a> {
    fn new(chunk: &'a [u8]) -> Self {
        Self { pos: 0, chunk }
    }
}

impl<'a> Iterator for SliceKoutputChunk<'a> {
    type Item = &'a [u8];
    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.chunk.len() {
            return None;
        }
        let out;
        if let Some(pos) = memchr(b'\n', &self.chunk[self.pos ..]) {
            out = &self.chunk[self.pos .. self.pos + pos];
            self.pos += pos + 1;
        } else {
            out = &self.chunk[self.pos ..];
            self.pos = self.chunk.len();
        }
        Some(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_lines(chunk: SliceKoutputChunk) -> Vec<&[u8]> {
        chunk.collect()
    }

    #[test]
    fn test_multiple_chunks_with_newlines() {
        let data = b"line1\nline2\nline3\nline4\nline5\n";
        let mut reader = SliceChunkReader::with_capacity(10, data);
        let mut results = Vec::new();

        while let Some(chunk) = reader.next() {
            results.extend(to_lines(chunk));
        }

        assert_eq!(results, vec![
            b"line1", b"line2", b"line3", b"line4", b"line5"
        ]);
    }

    #[test]
    fn test_single_chunk_without_trailing_newline() {
        let data = b"line1\nline2\nline3";
        let mut reader = SliceChunkReader::with_capacity(10, data);
        let mut results = Vec::new();

        while let Some(chunk) = reader.next() {
            results.extend(to_lines(chunk));
        }

        assert_eq!(results, vec![b"line1", b"line2", b"line3"]);
    }

    #[test]
    fn test_empty_input() {
        let data = b"";
        let mut reader = SliceChunkReader::with_capacity(10, data);
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_long_line_crossing_chunks() {
        let data = b"short1\naveryveryveryveryveryverylonglinewithoutnewline";
        let mut reader = SliceChunkReader::with_capacity(10, data);
        let mut results = Vec::new();

        while let Some(chunk) = reader.next() {
            results.extend(to_lines(chunk));
        }

        assert_eq!(results, vec![
            &b"short1"[..],
            &b"averyveryveryveryveryverylonglinewithoutnewline"[..]
        ]);
    }

    #[test]
    fn test_exact_chunk_size_split() {
        let data = b"abc\ndef\nghi\n";
        let mut reader = SliceChunkReader::with_capacity(8, data);
        let mut results = Vec::new();

        while let Some(chunk) = reader.next() {
            results.extend(to_lines(chunk));
        }

        assert_eq!(results, vec![b"abc", b"def", b"ghi"]);
    }
}
