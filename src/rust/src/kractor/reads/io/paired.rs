use std::fs::File;
use std::io::{BufWriter, Read};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use rustc_hash::FxHashSet as HashSet;

use crate::batchsender::BatchSender;
use crate::parser::fasta::FastaRecord;
use crate::parser::fastq::{FastqBytesChunkPairedReader, FastqContainer};
use crate::reader::bytes::BytesReader;

pub fn reader_kractor_paired_read<R1: Read + Send, R2: Read + Send>(
    id_sets: HashSet<&[u8]>,
    reader1: R1,
    ofile1: &str,
    reader2: R2,
    ofile2: &str,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Create output file and wrap in buffered writer
    let mut writer1 = BufWriter::with_capacity(buffer_size, File::create(ofile1)?);
    let mut writer2 = BufWriter::with_capacity(buffer_size, File::create(ofile2)?);

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<(FastaRecord<Bytes>, FastaRecord<Bytes>)>>,
            Receiver<Vec<(FastaRecord<Bytes>, FastaRecord<Bytes>)>>,
        ) = crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for (record1, record2) in chunk {
                    record1.write(&mut writer1)?;
                    record2.write(&mut writer2)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let mut reader1 = BytesReader::new(reader1);
        reader1.set_label("fq1");
        let mut reader2 = BytesReader::new(reader2);
        reader2.set_label("fq2");
        let mut reader = FastqBytesChunkPairedReader::with_capacity(chunk_size, reader1, reader2);

        let parser_handle = scope.spawn(move || {
            // will move `reader`, `parser_tx`, and `id_sets`
            // Shared atomic flag to signal if any thread encountered an error
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
            // Reuse the shared ID set for matching records
            let id_sets = &id_sets;
            // Rayon scope for spawning parsing threads
            rayon::scope(|s| -> Result<()> {
                while let Some(chunk) = reader.chunk_reader()? {
                    // If an error already occurred, stop spawning new threads
                    if has_error.load(Relaxed) {
                        return Ok(());
                    }
                    // Initialize a thread-local batch sender for matching records
                    let mut thread_tx = BatchSender::with_capacity(batch_size, parser_tx.clone());
                    // Clone the shared error state for this thread
                    let thread_has_error = has_error.clone();
                    let thread_err_tx = err_tx.clone();
                    s.spawn(move |_| {
                        let mut container1 = FastqContainer::new();
                        let mut container2 = FastqContainer::new();
                        loop {
                            match chunk.read_record(&mut container1, &mut container2) {
                                Ok(value) => match value {
                                    Some((record1, record2)) => {
                                        if id_sets.contains(&record1.id[..]) {
                                            // Attempt to send the matching record to the writer thread.
                                            match thread_tx.send((record1.into(), record2.into())) {
                                                Ok(_) => continue,
                                                Err(_) => {
                                                    // If this fails, it means the writer thread has already exited due to an error.
                                                    // Since that error will be reported separately, we can safely ignore this send failure.
                                                    // Mark that an error has occurred
                                                    thread_has_error.store(true, Relaxed);
                                                    return ();
                                                }
                                            };
                                        }
                                    }
                                    None => break,
                                },
                                Err(e) => {
                                    // Mark that an error has occurred
                                    thread_has_error.store(true, Relaxed);
                                    // Only the first error is reported
                                    let _ = thread_err_tx.try_send(e);
                                    return ();
                                }
                            }
                        }
                        let _ = thread_tx.flush();
                    });
                }
                Ok(())
            })?;
            // Clean up: close the error channel and drop parser sender
            drop(err_tx);
            drop(parser_tx);

            // Report the first error if any thread encountered one
            if has_error.load(Relaxed) {
                let err = err_rx.recv()?; // Safe unwrap because has_error is true
                Err(anyhow!(err))
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

#[cfg(test)]
mod tests {
    use std::fs::{self, File};
    use std::io::Write;

    use tempfile::tempdir;

    use super::*;

    fn write_temp_fastq(file_path: &str, records: &[(&str, &str)]) {
        let mut file = File::create(file_path).unwrap();
        for (id, seq) in records {
            writeln!(file, "@{}", id).unwrap();
            writeln!(file, "{}", seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", "I".repeat(seq.len())).unwrap();
        }
    }

    #[test]
    fn test_reader_kractor_paired_read() {
        // Create a temporary directory
        let dir = tempdir().unwrap();
        let fq1_path = dir.path().join("test_R1.fq");
        let fq2_path = dir.path().join("test_R2.fq");
        let ofile1_path = dir.path().join("out_R1.fq");
        let ofile2_path = dir.path().join("out_R2.fq");

        // Create paired records
        let records = vec![
            ("read1", "ACTGACTGACTG"),
            ("read2", "TGCATGCATGCA"),
            ("read3", "GGGGCCCCAAAA"),
        ];

        // Write input FASTQ files
        write_temp_fastq(fq1_path.to_str().unwrap(), &records);
        write_temp_fastq(fq2_path.to_str().unwrap(), &records);

        // Create an ID set containing only read2
        let id_set: HashSet<&[u8]> = vec![b"read2"].into_iter().map(|s| s.as_ref()).collect();

        // Run the paired reader
        reader_kractor_paired_read(
            id_set,
            File::open(fq1_path.to_str().unwrap()).unwrap(),
            ofile1_path.to_str().unwrap(),
            File::open(fq2_path.to_str().unwrap()).unwrap(),
            ofile2_path.to_str().unwrap(),
            2,
            4096,
            2,
            Some(10),
        )
        .unwrap();

        // Read and verify output
        let out_r1 = fs::read_to_string(ofile1_path).unwrap();
        let out_r2 = fs::read_to_string(ofile2_path).unwrap();
        assert!(out_r1.contains("read2"));
        assert!(out_r2.contains("read2"));
        assert!(!out_r1.contains("read1"));
        assert!(!out_r2.contains("read1"));
        assert!(!out_r1.contains("read3"));
        assert!(!out_r2.contains("read3"));
    }
}
