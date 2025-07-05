use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

use aho_corasick::AhoCorasick;
use anyhow::{anyhow, Result};
use memchr::{memrchr, Memchr};

use super::kractor_match_aho;
use crate::batchsender::BatchSender;
use crate::reader::slice::{SliceLineReader, SliceProgressBarReader};

pub fn mmap_kractor_koutput(
    include_aho: AhoCorasick,
    exclude_aho: Option<AhoCorasick>,
    reader: SliceProgressBarReader,
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
        let mut reader = KoutputSliceChunkReader::with_capacity(chunk_size, reader);
        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `matcher`
            let has_error = AtomicBool::new(false);
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            rayon::scope(|s| -> Result<()> {
                while let Some(chunk) = reader.chunk_reader() {
                    if has_error.load(Relaxed) {
                        return Ok(());
                    }
                    s.spawn(|_| {
                        let mut thread_tx =
                            BatchSender::with_capacity(batch_size, parser_tx.clone());
                        for line in chunk {
                            if kractor_match_aho(&include_aho, &exclude_aho, line) {
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

pub struct KoutputSliceChunkReader<'a> {
    chunk_size: usize,
    reader: SliceProgressBarReader<'a>,
}

impl<'a> KoutputSliceChunkReader<'a> {
    #[allow(dead_code)]
    pub fn new(reader: SliceProgressBarReader<'a>) -> Self {
        Self::with_capacity(8 * 1024, reader)
    }

    pub fn with_capacity(capacity: usize, reader: SliceProgressBarReader<'a>) -> Self {
        Self {
            chunk_size: capacity,
            reader,
        }
    }

    pub fn chunk_reader(&mut self) -> Option<SliceLineReader<'a, Memchr<'a>>> {
        // If we've reached or passed the end of the input buffer, stop reading
        let buf = match self.reader.take_slice(self.chunk_size) {
            Some(data) => data,
            None => return None,
        };
        let chunk;
        if buf.eof() {
            chunk = buf.as_slice();
        } else {
            let bytes = buf.as_slice();
            if let Some(pos) = memrchr(b'\n', bytes) {
                chunk = &bytes[..= pos];
            } else {
                self.chunk_size *= 2;
                return self.chunk_reader();
            }
        }
        self.reader.advance(chunk.len());
        Some(SliceLineReader::new(chunk))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn to_lines<'a>(mut chunk: SliceLineReader<'a, Memchr<'a>>) -> Vec<&'a [u8]> {
        let mut out = Vec::new();
        while let Some(slice) = chunk.read_line_inclusive() {
            out.push(slice);
        }
        out
    }

    #[test]
    fn test_multiple_chunks_with_newlines() {
        let data = b"line1\nline2\nline3\nline4\nline5\n";
        let mut reader =
            KoutputSliceChunkReader::with_capacity(10, SliceProgressBarReader::new(data));
        let mut results = Vec::new();

        while let Some(chunk) = reader.chunk_reader() {
            results.extend(to_lines(chunk));
        }

        assert_eq!(results, vec![
            b"line1\n", b"line2\n", b"line3\n", b"line4\n", b"line5\n"
        ]);
    }

    #[test]
    fn test_single_chunk_without_trailing_newline() {
        let data = b"line1\nline2\nline3";
        let mut reader =
            KoutputSliceChunkReader::with_capacity(10, SliceProgressBarReader::new(data));
        let mut results = Vec::new();

        while let Some(chunk) = reader.chunk_reader() {
            results.extend(to_lines(chunk));
        }

        assert_eq!(results, vec![
            &b"line1\n"[..],
            &b"line2\n"[..],
            &b"line3"[..]
        ]);
    }

    #[test]
    fn test_empty_input() {
        let data = b"";
        let mut reader =
            KoutputSliceChunkReader::with_capacity(10, SliceProgressBarReader::new(data));
        assert!(reader.chunk_reader().is_none());
    }

    #[test]
    fn test_long_line_crossing_chunks() {
        let data = b"short1\naveryveryveryveryveryverylonglinewithoutnewline";
        let mut reader =
            KoutputSliceChunkReader::with_capacity(10, SliceProgressBarReader::new(data));
        let mut results = Vec::new();

        while let Some(chunk) = reader.chunk_reader() {
            results.extend(to_lines(chunk));
        }

        assert_eq!(results, vec![
            &b"short1\n"[..],
            &b"averyveryveryveryveryverylonglinewithoutnewline"[..]
        ]);
    }

    #[test]
    fn test_exact_chunk_size_split() {
        let data = b"abc\ndef\nghi\n";
        let mut reader =
            KoutputSliceChunkReader::with_capacity(8, SliceProgressBarReader::new(data));
        let mut results = Vec::new();

        while let Some(chunk) = reader.chunk_reader() {
            results.extend(to_lines(chunk));
        }
        assert_eq!(results, vec![b"abc\n", b"def\n", b"ghi\n"]);
    }
}
