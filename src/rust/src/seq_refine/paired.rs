use std::io::{BufReader, Write};
use std::iter::zip;
use std::path::Path;

use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::{CompressionLvl, Compressor};

use crate::batchsender::BatchSender;
use crate::{fastq_reader::*, gz_compressed};
use crate::parser::fastq::FastqRecord;
use crate::{fastq_reader::FastqReader, seq_action::*};

pub(crate) fn reader_seq_refine_paired_read<P: AsRef<Path> + ?Sized>(
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
    chunk_size: usize,
    buffer_size: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let compression_level = CompressionLvl::new(compression_level)
        .map_err(|e| anyhow!("Invalid 'compression_level': {:?}", e))?;
    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (writer_tx, writer_rx): (Sender<(Option<Vec<u8>>, Option<Vec<u8>>)>, Receiver<(Option<Vec<u8>>, Option<Vec<u8>>)>) =
            crate::new_channel(nqueue);

        let (writer1_tx, writer1_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) =
            crate::new_channel(nqueue);

        let (writer2_tx, writer2_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) =
            crate::new_channel(nqueue);

        let mut reader1_tx_vec = Vec::with_capacity(threads);
        let mut reader1_rx_vec = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let (reader_tx, reader_rx): (
                Sender<Vec<FastqRecord<Bytes>>>,
                Receiver<Vec<FastqRecord<Bytes>>>,
            ) = crate::new_channel(nqueue);
            reader1_tx_vec.push(reader_tx);
            reader1_rx_vec.push(reader_rx);
        }
        let mut reader2_tx_vec = Vec::with_capacity(threads);
        let mut reader2_rx_vec = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let (reader_tx, reader_rx): (
                Sender<Vec<FastqRecord<Bytes>>>,
                Receiver<Vec<FastqRecord<Bytes>>>,
            ) = crate::new_channel(nqueue);
            reader2_tx_vec.push(reader_tx);
            reader2_rx_vec.push(reader_rx);
        }

        let (writer1_handle, gzip1) = if let Some(output_path) = output1_path {
            let output: &Path = output_path.as_ref();
            let handle = Some(scope.spawn(move || -> Result<()> {
                let mut writer = fastq_writer(output, output1_bar)?;
                for chunk in writer1_rx {
                    writer.write_all(&chunk).map_err(|e| {
                        anyhow!("(Writer1) Failed to write FastqRecord to output: {}", e)
                    })?;
                }
                writer
                    .flush()
                    .map_err(|e| anyhow!("(Writer1) Failed to flush writer:  {}", e))?;
                Ok(())
            }));
            let gzip = gz_compressed(output);
            (handle, gzip)
        } else {
            (None, false)
        };

        let (writer2_handle, gzip2) = if let Some(output_path) = output2_path {
            let output: &Path = output_path.as_ref();
            let handle =Some(scope.spawn(move || -> Result<()> {
                let mut writer = fastq_writer(output, output2_bar)?;
                for chunk in writer2_rx {
                    writer.write_all(&chunk).map_err(|e| {
                        anyhow!("(Writer2) Failed to write FastqRecord to output: {}", e)
                    })?;
                }
                writer
                    .flush()
                    .map_err(|e| anyhow!("(Writer2) Failed to flush writer: {}", e))?;
                Ok(())
            }));
            let gzip = gz_compressed(output);
            (handle, gzip)
        } else {
            (None, false)
        };

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for (records1, records2) in writer_rx {
                if let Some(records1) = records1 {
                    writer1_tx.send(records1).map_err(|e| {
                        anyhow!(
                            "(Writer dispatch) Failed to send read1 batch to Writer1 thread: {}",
                            e
                        )
                    })?;
                }
                if let Some(records2) = records2 {
                    writer2_tx.send(records2).map_err(|e| {
                        anyhow!(
                            "(Writer dispatch) Failed to send read2 batch to Writer2 thread: {}",
                            e
                        )
                    })?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let has_writer1 = writer2_handle.is_some();
        let has_writer2 = writer2_handle.is_some();
        let mut parser_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let reader1_rx = reader1_rx_vec.pop().ok_or(anyhow!(
                "(Parser setup) Not enough read1 channels for threads"
            ))?;
            let reader2_rx = reader2_rx_vec.pop().ok_or(anyhow!(
                "(Parser setup) Not enough read2 channels for threads"
            ))?;
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                let mut records1_pool: Vec<u8> = Vec::with_capacity(crate::BLOCK_SIZE);
                let mut records2_pool: Vec<u8> = Vec::with_capacity(crate::BLOCK_SIZE);
                let mut compressor = Compressor::new(compression_level);
                loop {
                    let (records1, records2) = match (reader1_rx.recv(), reader2_rx.recv()) {
                        (Ok(rec1), Ok(rec2)) => (rec1, rec2),
                        (Err(_), Ok(_)) => {
                            return Err(anyhow!(
                                "(Parser) FASTQ pairing error: read1 channel closed before read2"
                            ));
                        }
                        (Ok(_), Err(_)) => {
                            return Err(anyhow!(
                                "(Parser) FASTQ pairing error: read2 channel closed before read1"
                            ));
                        }
                        (Err(_), Err(_)) => {
                            break;
                        }
                    };
                    if records1.len() != records2.len() {
                        return Err(anyhow!("(Parser) FASTQ pairing error: record count mismatch (read1: {}, read2: {})", records1.len(), records2.len()));
                    }

                    // Initialize a thread-local batch sender for matching records
                    for (mut record1, mut record2) in zip(records1, records2) {
                        actions.transform_fastq_bytes(&mut record1, &mut record2)?;
                        if (records1_pool.capacity() - records1_pool.len() < record1.bytes_size()) ||
                            (records2_pool.capacity() - records2_pool.len() < record2.bytes_size()) {
                            let pack1 = if has_writer1 {
                                Some(fastq_pack(&records1_pool, &mut compressor, gzip1)?)
                            } else {
                                None
                            };
                            let pack2 = if has_writer2 {
                                Some(fastq_pack(&records2_pool, &mut compressor, gzip2)?)
                            } else {
                                None
                            };
                            tx.send((pack1, pack2)).map_err(|e| {
                                anyhow!(
                                    "(Parser) Failed to send send parsed record pair to Writer thread: {}",
                                    e
                                )
                            })?;
                            records1_pool.clear();
                            records2_pool.clear();
                        }
                        record1.extend(&mut records1_pool);
                        record2.extend(&mut records2_pool);
                    }
                }
                if !records1_pool.is_empty() {
                    let pack1 = if has_writer1 {
                        Some(fastq_pack(&records1_pool, &mut compressor, gzip1)?)
                    } else {
                        None
                    };
                    let pack2 = if has_writer2 {
                        Some(fastq_pack(&records2_pool, &mut compressor, gzip2)?)
                    } else {
                        None
                    };
                    tx.send((pack1, pack2)).map_err(|e| {
                        anyhow!(
                            "(Parser) Failed to send send parsed record pair to Writer thread: {}",
                            e
                        )
                    })?;
                }
                Ok(())
            });
            parser_handles.push(handle);
        }
        drop(writer_tx);

        // ─── reader Thread ─────────────────────────────────────
        let input1: &Path = input1_path.as_ref();
        let reader1_handle = scope.spawn(move || -> Result<()> {
            let reader = fastq_reader(input1, input1_bar)?;
            let mut reader = FastqReader::new(BufReader::with_capacity(buffer_size, reader));
            let mut reader1_tx_vec = reader1_tx_vec
                .into_iter()
                .map(|tx| BatchSender::with_capacity(chunk_size, tx))
                .collect::<Vec<_>>();
            let mut n = 0usize;
            while let Some(record) = reader
                .read_record()
                .map_err(|e| anyhow!("(Reader1) Error while reading FASTQ record: {}", e))?
            {
                let idx = n % threads;
                reader1_tx_vec[idx].send(record).map_err(|e| {
                    anyhow!(
                        "(Reader1) Failed to send FASTQ record to Parser thread: {}",
                        e
                    )
                })?;
                n += 1;
            }
            for mut reader_tx in reader1_tx_vec {
                reader_tx.flush().map_err(|e| {
                    anyhow!("(Reader1) Failed to flush records to Parser thread: {}", e)
                })?;
            }
            Ok(())
        });

        let input2: &Path = input2_path.as_ref();
        let reader2_handle = scope.spawn(move || -> Result<()> {
            let reader2 = fastq_reader(input2, input2_bar)?;
            let mut reader2 = FastqReader::new(BufReader::with_capacity(buffer_size, reader2));
            let mut reader2_tx_vec = reader2_tx_vec
                .into_iter()
                .map(|tx| BatchSender::with_capacity(chunk_size, tx))
                .collect::<Vec<_>>();
            let mut n = 0usize;
            while let Some(record) = reader2
                .read_record()
                .map_err(|e| anyhow!("(Reader2) Error while reading FASTQ record: {}", e))?
            {
                let idx = n % threads;
                reader2_tx_vec[idx].send(record).map_err(|e| {
                    anyhow!(
                        "(Reader2) Failed to send FASTQ record to Parser thread: {}",
                        e
                    )
                })?;
                n += 1;
            }
            for mut reader_tx in reader2_tx_vec {
                reader_tx.flush().map_err(|e| {
                    anyhow!("(Reader2) Failed to flush records to Parser thread: {}", e)
                })?;
            }
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
            .map_err(|e| anyhow!("(Writer) thread panicked: {:?}", e))??;

        for handler in parser_handles {
            handler
                .join()
                .map_err(|e| anyhow!("(Parser) thread panicked: {:?}", e))??;
        }
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
        let ranges: SortedSeqRanges = vec![SeqRange::From(3), SeqRange::To(2)]
            .into_iter()
            .collect();
        let mut actions = SubseqActions::builder();
        actions.add_action(SeqAction::Embed("UMI".to_string()), ranges);
        let paired_actions = SubseqPairedActions::new(Some(actions.build()), None);

        // Run paired reader pipeline
        reader_seq_refine_paired_read(
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
