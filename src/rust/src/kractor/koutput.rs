use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::sync::atomic::{AtomicBool, Ordering};

use aho_corasick::{AhoCorasick, AhoCorasickKind};
use anyhow::{anyhow, Result};
use bytes::{Bytes, BytesMut};
use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use memchr::{memchr, memmem, memrchr};
use memmap2::{Advice, Mmap};
use rustc_hash::FxHashSet as HashSet;

use crate::kractor::batchsender::BatchSender;

fn new_channel<T>(queue: Option<usize>) -> (Sender<T>, Receiver<T>) {
    if let Some(queue) = queue {
        bounded(queue)
    } else {
        unbounded()
    }
}

#[allow(clippy::too_many_arguments)]
pub fn reader_kractor_koutput(
    file: &str,
    ofile: &str,
    patterns: &[&str],
    batch_size: usize,
    read_buffer: usize,
    write_buffer: usize,
    read_queue: Option<usize>,
    write_queue: Option<usize>,
) -> Result<()> {
    let mut reader = std::fs::File::open(file)?;
    let mut writer =
        BufWriter::with_capacity(write_buffer, std::fs::File::create(ofile)?);
    let matcher = AhoCorasick::builder()
        .kind(Some(AhoCorasickKind::DFA))
        .build(patterns)?;

    std::thread::scope(|scope| {
        let (reader_tx, parser_rx): (Sender<Bytes>, Receiver<Bytes>) =
            new_channel(read_queue);
        let (parser_tx, writer_rx): (Sender<Vec<Bytes>>, Receiver<Vec<Bytes>>) =
            new_channel(write_queue);

        // first thread, used to write lines
        let writer_handle = scope.spawn(move || -> Result<()> {
            for chunk in writer_rx {
                for line in chunk {
                    writer.write_all(&line)?;
                }
            }
            Ok(())
        });

        // second thread, used to filter lines
        let parser_handle = scope.spawn(|| {
            let matcher_ref = &matcher;
            let parser_batch_tx =
                BatchSender::with_capacity(batch_size, parser_tx);
            rayon::scope(|s| -> Result<()> {
                let mut main_tx = parser_batch_tx.clone();
                let mut leftover: Option<Bytes> = None;
                let mut start;
                for chunk in parser_rx {
                    // merge the leftover with current chunk
                    if let Some(head) = leftover {
                        // Combine previous partial line with new chunk — 1 allocation, unavoidable
                        if let Some(pos) = memchr(b'\n', &chunk) {
                            let line: Bytes = head
                                .iter()
                                .chain(&chunk[..= pos])
                                .copied()
                                .collect();
                            if kractor_match_aho(matcher_ref, &line) {
                                main_tx.send(line)?;
                            };
                            start = pos + 1;
                        } else {
                            let line: Bytes = head
                                .iter()
                                .chain(chunk.iter())
                                .copied()
                                .collect();
                            leftover = Some(line);
                            continue;
                        }
                    } else {
                        start = 0;
                    }

                    // then we identify the chunk with multiple records-in koutput file,
                    // each line is a single record
                    if let Some(pos) = memrchr(b'\n', &chunk[start ..]) {
                        let records = chunk.slice(start ..= start + pos);
                        let mut tx = parser_batch_tx.clone();
                        s.spawn(move |_| {
                            let mut line_start = 0;
                            while let Some(pos) =
                                memchr(b'\n', &records[line_start ..])
                            {
                                let line = records
                                    .slice(line_start ..= line_start + pos);
                                if kractor_match_aho(matcher_ref, &line) {
                                    match tx.send(line) {
                                        // Only error when `writer_handle` return error
                                        Err(_) => {
                                            return ();
                                        }
                                        Ok(_) => {}
                                    };
                                }
                                line_start += pos + 1
                            }
                            let _ = tx.flush();
                        });
                        start += pos + 1;
                    }
                    // handle trailing bytes (no newline)
                    if start < chunk.len() {
                        leftover = Some(chunk.slice(start ..)); // save remainder (no copy)
                    } else {
                        leftover = None
                    }
                }
                if let Some(line) = leftover {
                    if kractor_match_aho(matcher_ref, &line) {
                        main_tx.send(Bytes::from(line))?;
                    }
                }
                main_tx.flush()?;
                Ok(())
            })
        });

        // third thread, used to read lines
        let reader_handle = scope.spawn(move || -> Result<()> {
            let mut buf;
            loop {
                buf = BytesMut::with_capacity(read_buffer);
                // read files and pass chunks to parser
                unsafe { buf.set_len(read_buffer) };
                let nbytes = reader.read(&mut buf)?;
                if nbytes == 0 {
                    break;
                }
                if nbytes < buf.len() {
                    // SAFETY: Shrinking the buffer cannot expose uninitialized bytes.
                    unsafe { buf.set_len(nbytes) };
                }
                reader_tx.send(buf.freeze())?;
            }
            Ok(())
        });
        reader_handle
            .join()
            .map_err(|e| anyhow!("Reader thread panicked: {:?}", e))??;
        let parser_result = parser_handle
            .join()
            .map_err(|e| anyhow!("Parser thread panicked: {:?}", e))?;
        writer_handle
            .join()
            .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))??;
        parser_result?;
        Ok(())
    })
}

pub fn mmap_kractor_koutput(
    file: &str,
    ofile: &str,
    patterns: &[&str],
    batch_size: usize,
    write_buffer: usize,
    read_queue: Option<usize>,
    write_queue: Option<usize>,
) -> Result<()> {
    let file = File::open(file)?;
    let matcher = AhoCorasick::builder()
        .kind(Some(AhoCorasickKind::DFA))
        .build(patterns)?;

    // let taxid_sets: HashSet<&[u8]> =
    //     patterns.into_iter().map(|s| s.as_bytes()).collect();
    // https://github.com/rayon-rs/rayon/discussions/1164
    let map = unsafe { Mmap::map(&file) }?;
    map.advise(Advice::Sequential)?;

    let mut writer =
        BufWriter::with_capacity(write_buffer, std::fs::File::create(ofile)?);
    std::thread::scope(|scope| {
        let (reader_tx, parser_rx) = new_channel(read_queue);
        let (parser_tx, writer_rx) = new_channel(write_queue);

        // first thread, used to write lines
        let writer_handle = scope.spawn(|| -> Result<()> {
            for chunk in writer_rx {
                for line in chunk {
                    writer.write_all(line)?;
                }
            }
            Ok(())
        });

        // second thread, used to filter lines
        let parser_handle = scope.spawn(|| {
            let matcher_ref = &matcher;
            let parser_batch_tx =
                BatchSender::with_capacity(batch_size, parser_tx);
            rayon::scope(|s| {
                for chunk in parser_rx {
                    // let taxid_sets_ref = &taxid_sets;
                    let mut tx = parser_batch_tx.clone();
                    s.spawn(move |_| {
                        for line in chunk {
                            if kractor_match_aho(matcher_ref, line) {
                                match tx.send(line) {
                                    // Only error when `writer_handle` return error
                                    // So we just omit the error information
                                    Err(_) => {
                                        return ();
                                    }
                                    Ok(_) => continue,
                                };
                            }
                        }
                        let _ = tx.flush();
                    });
                }
            });
        });

        // third thread, used to read lines
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
            // the last line may not contain `\n`
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
            .map_err(|e| anyhow!("Parser thread panicked: {:?}", e))?;
        writer_handle
            .join()
            .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))??;
        Ok(())
    })
}

static KOUTPUT_TAXID_PREFIX_FINDER: std::sync::LazyLock<memmem::Finder> =
    std::sync::LazyLock::new(|| memmem::Finder::new("(taxid"));
static KOUTPUT_TAXID_PREFIX_FINDERREV: std::sync::LazyLock<memmem::FinderRev> =
    std::sync::LazyLock::new(|| memmem::FinderRev::new("(taxid"));

#[allow(dead_code)]
fn kractor_match_aho(matcher: &AhoCorasick, line: &[u8]) -> bool {
    // println!("Matching line: {:?}", String::from_utf8_lossy(line));
    // Efficient 3rd column parsing
    let mut field_start = 0usize;
    let mut field_index = 0usize;
    while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
        if field_index == 2 {
            // we don't include the last `\t`
            let taxid = &line[field_start .. (field_start + tab_pos)];
            if let Some(start) = KOUTPUT_TAXID_PREFIX_FINDER.find(taxid) {
                let mut input = aho_corasick::Input::new(taxid);
                input.set_start(start);
                return matcher.find(taxid).is_some();
            } else {
                return false;
            }
        }
        field_index += 1;
        field_start += tab_pos + 1;
    }
    false
}

#[allow(dead_code)]
fn kractor_match_hash(taxid_sets: &HashSet<&[u8]>, line: &[u8]) -> bool {
    // Efficient 3rd column parsing
    let mut field_start = 0usize;
    let mut field_index = 0usize;
    while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
        if field_index == 2 {
            // we don't include the last `\t`
            let taxa = &line[field_start .. (field_start + tab_pos)];
            // extract the taxid pattern
            if let Some(start) = KOUTPUT_TAXID_PREFIX_FINDERREV.rfind(taxa) {
                if let Some(pos) = memchr(b')', &taxa[start ..]) {
                    return taxid_sets.contains(&taxa[start ..= start + pos]);
                };
            };
            return false;
        }
        field_start += tab_pos + 1;
        field_index += 1;
    }
    false
}

#[cfg(test)]
mod tests {
    use std::fs::{read_to_string, File};
    use std::io::Write;

    use tempfile::tempdir;

    use super::*;

    fn write_sample_koutput(path: &std::path::Path) -> Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "read1\tumi1\t(taxid12345)\tother")?;
        writeln!(file, "read2\tumi2\t(taxid54321)\tother")?;
        writeln!(file, "read3\tumi3\t(no_taxid)\tother")?;
        writeln!(file, "read4\tumi4\t(taxid99999)\tother")?;
        Ok(())
    }

    #[test]
    fn test_reader_kractor_koutput() -> Result<()> {
        let dir = tempdir()?;
        let input_path = dir.path().join("input.txt");
        let output_path = dir.path().join("output.txt");

        write_sample_koutput(&input_path)?;

        let patterns = &["taxid12345", "taxid99999"];
        reader_kractor_koutput(
            input_path.to_str().unwrap(),
            output_path.to_str().unwrap(),
            patterns,
            2,
            64,
            64,
            Some(10),
            Some(10),
        )?;

        let output = read_to_string(output_path)?;
        assert!(output.contains("taxid12345"));
        assert!(output.contains("taxid99999"));
        assert!(!output.contains("taxid54321"));
        assert!(!output.contains("no_taxid"));
        Ok(())
    }

    #[test]
    fn test_mmap_kractor_koutput() -> Result<()> {
        let dir = tempdir()?;
        let input_path = dir.path().join("input.txt");
        let output_path = dir.path().join("output.txt");

        write_sample_koutput(&input_path)?;

        let patterns = &["taxid54321"];
        mmap_kractor_koutput(
            input_path.to_str().unwrap(),
            output_path.to_str().unwrap(),
            patterns,
            2,
            64,
            Some(10),
            Some(10),
        )?;

        let output = read_to_string(output_path)?;
        assert!(output.contains("taxid54321"));
        assert!(!output.contains("taxid12345"));
        assert!(!output.contains("taxid99999"));
        Ok(())
    }
}
