use std::io::BufWriter;
use std::io::Write;
use std::iter::zip;
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::{CompressionLvl, Compressor};

use crate::batchsender::BatchSender;
use crate::fastq_reader::*;
use crate::parser::fastq::FastqRecord;
use crate::seq_action::*;
use crate::utils::*;

pub(crate) fn seq_refine_paired_read<P: AsRef<Path> + ?Sized>(
    input1_path: &P,
    input1_bar: Option<ProgressBar>,
    input2_path: &P,
    input2_bar: Option<ProgressBar>,
    output1_path: Option<&P>,
    output1_bar: Option<ProgressBar>,
    output2_path: Option<&P>,
    output2_bar: Option<ProgressBar>,
    actions: &SubseqPairedActions,
    compression_level: i32,
    batch_size: usize,
    chunk_bytes: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let compression_level = CompressionLvl::new(compression_level)
        .map_err(|e| anyhow!("Invalid 'compression_level': {:?}", e))?;
    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (writer_tx, writer_rx): (
            Sender<(Option<Vec<u8>>, Option<Vec<u8>>)>,
            Receiver<(Option<Vec<u8>>, Option<Vec<u8>>)>,
        ) = new_channel(nqueue);
        let (writer1_tx, writer1_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = new_channel(nqueue);
        let (writer2_tx, writer2_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = new_channel(nqueue);

        let (reader_tx, reader_rx): (
            Sender<(Vec<FastqRecord<Bytes>>, Vec<FastqRecord<Bytes>>)>,
            Receiver<(Vec<FastqRecord<Bytes>>, Vec<FastqRecord<Bytes>>)>,
        ) = new_channel(nqueue);
        let (reader1_tx, reader1_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = new_channel(nqueue);
        let (reader2_tx, reader2_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        let (writer1_handle, gzip1) = if let Some(output_path) = output1_path {
            let output: &Path = output_path.as_ref();
            let handle = Some(scope.spawn(move || -> Result<()> {
                let mut writer =
                    BufWriter::with_capacity(chunk_bytes, new_writer(output, output1_bar)?);
                for chunk in writer1_rx {
                    writer.write_all(&chunk).with_context(|| {
                        format!("(Writer1) Failed to write Fastq records to output")
                    })?;
                }
                writer
                    .flush()
                    .with_context(|| format!("(Writer1) Failed to flush writer"))?;
                Ok(())
            }));
            let gzip = gz_compressed(output);
            (handle, gzip)
        } else {
            (None, false)
        };

        let (writer2_handle, gzip2) = if let Some(output_path) = output2_path {
            let output: &Path = output_path.as_ref();
            let handle = Some(scope.spawn(move || -> Result<()> {
                let mut writer =
                    BufWriter::with_capacity(chunk_bytes, new_writer(output, output2_bar)?);
                for chunk in writer2_rx {
                    writer.write_all(&chunk).with_context(|| {
                        format!("(Writer2) Failed to write Fastq records to output")
                    })?;
                }
                writer
                    .flush()
                    .with_context(|| format!("(Writer2) Failed to flush writer"))?;
                Ok(())
            }));
            let gzip = gz_compressed(output);
            (handle, gzip)
        } else {
            (None, false)
        };

        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for (records1, records2) in writer_rx {
                if let Some(records1) = records1 {
                    writer1_tx.send(records1).with_context(|| {
                        format!("(Writer dispatch) Failed to send read1 batch to Writer1 thread")
                    })?;
                }
                if let Some(records2) = records2 {
                    writer2_tx.send(records2).with_context(|| {
                        format!("(Writer dispatch) Failed to send read2 batch to Writer2 thread")
                    })?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let has_writer1 = writer1_handle.is_some();
        let has_writer2 = writer2_handle.is_some();
        let mut parser_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                let mut records1_pool: Vec<u8> = Vec::with_capacity(chunk_bytes);
                let mut records2_pool: Vec<u8> = Vec::with_capacity(chunk_bytes);
                let mut compressor = Compressor::new(compression_level);
                while let Ok((records1, records2)) = rx.recv() {
                    // Initialize a thread-local batch sender for matching records
                    for (mut record1, mut record2) in zip(records1, records2) {
                        if record1.id != record2.id {
                            return Err(anyhow!(
                                "FASTQ pairing error: record1 ID = {}, record2 ID = {}. These records do not match and cannot be paired.",
                                String::from_utf8_lossy(&record1.id),
                                String::from_utf8_lossy(&record2.id)
                            ));
                        }
                        actions.transform_fastq(&mut record1, &mut record2)?;
                        if records1_pool.capacity() - records1_pool.len() < record1.bytes_size() ||
                            records2_pool.capacity() - records2_pool.len() < record2.bytes_size() {
                            let pack1 = if has_writer1 {
                                let mut pack = Vec::with_capacity(chunk_bytes);
                                std::mem::swap(&mut records1_pool, &mut pack);
                                if gzip1 {
                                    pack = gzip_pack(&pack, &mut compressor)?
                                }
                                Some(pack)
                            } else {
                                None
                            };
                            let pack2 = if has_writer2 {
                                let mut pack = Vec::with_capacity(chunk_bytes);
                                std::mem::swap(&mut records2_pool, &mut pack);
                                if gzip2 {
                                    pack = gzip_pack(&pack, &mut compressor)?
                                }
                                Some(pack)
                            } else {
                                None
                            };
                            tx.send((pack1, pack2)).with_context(|| {
                                format!(
                                    "(Parser) Failed to send send parsed record pair to Writer thread"
                                )
                            })?;
                        }
                        record1.extend(&mut records1_pool);
                        record2.extend(&mut records2_pool);
                    }
                }
                if !records1_pool.is_empty() {
                    let pack1 = if has_writer1 {
                        let pack = if gzip1 {
                            gzip_pack(&records1_pool, &mut compressor)?
                        } else {
                            records1_pool
                        };
                        Some(pack)
                    } else {
                        None
                    };
                    let pack2 = if has_writer2 {
                        let pack = if gzip2 {
                            gzip_pack(&records2_pool, &mut compressor)?
                        } else {
                            records2_pool
                        };
                        Some(pack)
                    } else {
                        None
                    };
                    tx.send((pack1, pack2)).with_context(|| {
                        format!(
                            "(Parser) Failed to send send parsed record pair to Writer thread"
                        )
                    })?;
                }
                Ok(())
            });
            parser_handles.push(handle);
        }
        drop(reader_rx);
        drop(writer_tx);

        // ─── reader Thread ─────────────────────────────────────
        let reader_handle = scope.spawn(move || -> Result<()> {
            loop {
                let (records1, records2) = match (reader1_rx.recv(), reader2_rx.recv()) {
                    (Ok(rec1), Ok(rec2)) => (rec1, rec2),
                    (Err(_), Ok(_)) => {
                        return Err(anyhow!(
                            "(Reader collect) FASTQ pairing error: read1 channel closed before read2"
                        ));
                    }
                    (Ok(_), Err(_)) => {
                        return Err(anyhow!(
                            "(Reader collect) FASTQ pairing error: read2 channel closed before read1"
                        ));
                    }
                    (Err(_), Err(_)) => {
                        break;
                    }
                };
                if records1.len() != records2.len() {
                    return Err(anyhow!("(Reader collect) FASTQ pairing error: record count mismatch (read1: {}, read2: {})", records1.len(), records2.len()));
                }
                reader_tx.send((records1, records2)).with_context(|| {
                    format!(
                        "(Reader collect) Failed to send send parsed record pair to Parser thread"
                    )
                })?;
            }
            Ok(())
        });

        let input1: &Path = input1_path.as_ref();
        let reader1_handle = scope.spawn(move || -> Result<()> {
            let mut reader = FastqReader::with_capacity(
                BUFFER_SIZE,
                new_reader(input1, BUFFER_SIZE, input1_bar)?,
            );
            let mut thread_tx = BatchSender::with_capacity(batch_size, reader1_tx);
            while let Some(record) = reader
                .read_record()
                .with_context(|| format!("(Reader1) Failed to read FASTQ record"))?
            {
                thread_tx.send(record).with_context(|| {
                    format!("(Reader1) Failed to send FASTQ record to reader collect thread")
                })?;
            }
            thread_tx.flush().with_context(|| {
                format!("(Reader1) Failed to flush records to reader collect thread")
            })?;
            Ok(())
        });

        let input2: &Path = input2_path.as_ref();
        let reader2_handle = scope.spawn(move || -> Result<()> {
            let mut reader = FastqReader::with_capacity(
                BUFFER_SIZE,
                new_reader(input2, BUFFER_SIZE, input2_bar)?,
            );
            let mut thread_tx = BatchSender::with_capacity(batch_size, reader2_tx);
            while let Some(record) = reader
                .read_record()
                .with_context(|| format!("(Reader2) Failed to read FASTQ record"))?
            {
                thread_tx.send(record).with_context(|| {
                    format!("(Reader2) Failed to send FASTQ record to reader collect thread")
                })?;
            }
            thread_tx.flush().with_context(|| {
                format!("(Reader2) Failed to flush records to reader collect thread")
            })?;
            Ok(())
        });

        // ─── Join Threads and Propagate Errors ────────────────
        if let Some(writer_handle) = writer1_handle {
            writer_handle
                .join()
                .map_err(|e| anyhow!("(Writer1) thread panicked: {:?}", e))??;
        };
        if let Some(writer_handle) = writer2_handle {
            writer_handle
                .join()
                .map_err(|e| anyhow!("(Writer2) thread panicked: {:?}", e))??;
        };
        writer_handle
            .join()
            .map_err(|e| anyhow!("(Writer dispatch) thread panicked: {:?}", e))??;

        for handler in parser_handles {
            handler
                .join()
                .map_err(|e| anyhow!("(Parser) thread panicked: {:?}", e))??;
        }
        reader_handle
            .join()
            .map_err(|e| anyhow!("(Reader collect) thread panicked: {:?}", e))??;
        reader1_handle
            .join()
            .map_err(|e| anyhow!("(Reader1) thread panicked: {:?}", e))??;
        reader2_handle
            .join()
            .map_err(|e| anyhow!("(Reader2) thread panicked: {:?}", e))??;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::Read;

    use isal::read::GzipDecoder;

    use super::*;
    use crate::seq_action::{SeqAction, SubseqActions};
    use crate::seq_range::{SeqRange, SeqRanges};

    #[test]
    fn test_reader_seq_refine_paired_read_passthrough() -> Result<()> {
        use std::fs::write;

        // Sample paired-end FASTQ records (matching IDs)
        let read1 = b"@SEQ_ID1\nACGT\n+\n!!!!\n@SEQ_ID2\nTGCA\n+\n####\n";
        let read2 = b"@SEQ_ID1\nTTAA\n+\n$$$$\n@SEQ_ID2\nAATT\n+\n%%%%\n";

        // Create temp directory and files
        let tmp = tempfile::tempdir()?;
        let in1_path = tmp.path().join("read1.fq");
        let in2_path = tmp.path().join("read2.fq");
        let out1_path = tmp.path().join("out1.fq.gz");
        let out2_path = tmp.path().join("out2.fq.gz");

        write(&in1_path, read1)?;
        write(&in2_path, read2)?;

        // UMI extraction rule
        let ranges: SeqRanges = vec![SeqRange::From(3), SeqRange::To(2)]
            .into_iter()
            .collect();
        let mut actions = SubseqActions::builder();
        actions
            .add_action(
                SeqAction::Embed(Bytes::from_owner("UMI".as_bytes())),
                ranges,
            )
            .unwrap();
        let paired_actions = SubseqPairedActions::new(Some(actions.build().unwrap()), None);

        // Run paired reader pipeline
        seq_refine_paired_read(
            &in1_path,
            None,
            &in2_path,
            None,
            Some(&out1_path),
            None,
            Some(&out2_path),
            None,
            &paired_actions,
            4,         // compression
            1,         // chunk size
            64 * 1024, // buffer size
            Some(2),   // queue size
            1,         // threads
        )?;

        // Decompress and read output
        let mut buf1 = String::new();
        GzipDecoder::new(File::open(out1_path)?).read_to_string(&mut buf1)?;
        let mut buf2 = String::new();
        GzipDecoder::new(File::open(out2_path)?).read_to_string(&mut buf2)?;

        // Check contents
        assert_eq!(
            buf1.as_bytes(),
            b"@SEQ_ID1 RSAHMI{UMI:ACT}\nACGT\n+\n!!!!\n@SEQ_ID2 RSAHMI{UMI:TGA}\nTGCA\n+\n####\n"
        );
        assert_eq!(
            buf2.as_bytes(),
            b"@SEQ_ID1 RSAHMI{UMI:ACT}\nTTAA\n+\n$$$$\n@SEQ_ID2 RSAHMI{UMI:TGA}\nAATT\n+\n%%%%\n"
        );

        Ok(())
    }
}
