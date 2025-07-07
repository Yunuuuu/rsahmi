use std::io::{BufReader, Read, Write};
use std::iter::zip;

use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};

use crate::batchsender::BatchSender;
use crate::parser::fastq::FastqRecord;
use crate::{fastq_reader::FastqReader, seq_action::*};

pub(crate) fn reader_seq_refine_paired_read<
    R1: Read + Send,
    R2: Read + Send,
    W1: Write + Send,
    W2: Write + Send,
>(
    reader1: &mut R1,
    writer1: &mut Option<W1>,
    reader2: &mut R2,
    writer2: &mut Option<W2>,
    actions: &SubseqPairedActions,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let mut reader1 = FastqReader::new(BufReader::with_capacity(buffer_size, reader1));
    let mut reader2 = FastqReader::new(BufReader::with_capacity(buffer_size, reader2));
    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (writer_tx, writer_rx): (
            Sender<Vec<(FastqRecord<Bytes>, FastqRecord<Bytes>)>>,
            Receiver<Vec<(FastqRecord<Bytes>, FastqRecord<Bytes>)>>,
        ) = crate::new_channel(nqueue);

        let (writer1_tx, writer1_rx): (Sender<Vec<Vec<u8>>>, Receiver<Vec<Vec<u8>>>) =
            crate::new_channel(nqueue);

        let (writer2_tx, writer2_rx): (Sender<Vec<Vec<u8>>>, Receiver<Vec<Vec<u8>>>) =
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

        let has_writer1 = writer1.is_some();
        let writer1_handle = if let Some(writer) = writer1 {
            Some(scope.spawn(move || -> Result<()> {
                for chunk in writer1_rx {
                    for record in chunk {
                        writer.write_all(&record).map_err(|e| {
                            anyhow!("(Writer1) Failed to write FastqRecord to output: {}", e)
                        })?;
                    }
                }
                writer
                    .flush()
                    .map_err(|e| anyhow!("(Writer1) Failed to flush writer:  {}", e))?;
                Ok(())
            }))
        } else {
            None
        };

        let has_writer2 = writer2.is_some();
        let writer2_handle = if let Some(writer) = writer2 {
            Some(scope.spawn(move || -> Result<()> {
                for chunk in writer2_rx {
                    for record in chunk {
                        writer.write_all(&record).map_err(|e| {
                            anyhow!("(Writer2) Failed to write FastqRecord to output: {}", e)
                        })?;
                    }
                }
                writer
                    .flush()
                    .map_err(|e| anyhow!("(Writer2) Failed to flush writer: {}", e))?;
                Ok(())
            }))
        } else {
            None
        };

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                let (records1, records2): (Vec<FastqRecord<Bytes>>, Vec<FastqRecord<Bytes>>) =
                    chunk.into_iter().unzip();
                if has_writer1 {
                    let records1: Vec<Vec<u8>> = records1
                        .into_iter()
                        .map(|recorde| recorde.as_vec())
                        .collect();
                    writer1_tx.send(records1).map_err(|e| {
                        anyhow!(
                            "(Writer dispatch) Failed to send read1 batch to Writer1 thread: {}",
                            e
                        )
                    })?;
                }
                if has_writer2 {
                    let records2: Vec<Vec<u8>> = records2
                        .into_iter()
                        .map(|recorde| recorde.as_vec())
                        .collect();
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
                let mut thread_tx = BatchSender::with_capacity(batch_size, tx);
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
                        thread_tx
                            .send((record1, record2))
                            .map_err(|e| anyhow!("(Parser) Failed to send parsed record pair to Writer thread: {}", e))?;
                    }
                }
                thread_tx
                    .flush()
                    .map_err(|e| anyhow!("(Parser) Failed to flush remaining records to Writer thread: {}", e))?;
                Ok(())
            });
            parser_handles.push(handle);
        }
        drop(writer_tx);

        // ─── reader Thread ─────────────────────────────────────
        let reader1_handle = scope.spawn(move || -> Result<()> {
            let mut reader1_tx_vec = reader1_tx_vec
                .into_iter()
                .map(|tx| BatchSender::with_capacity(chunk_size, tx))
                .collect::<Vec<_>>();
            let mut n = 0usize;
            while let Some(record) = reader1
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

        let reader2_handle = scope.spawn(move || -> Result<()> {
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

    use flate2::read::GzDecoder;
    use tempfile::tempdir;

    use super::*;
    use crate::fastq_reader::fastq_writer;

    #[test]
    fn test_reader_seq_refine_paired_read_passthrough() -> Result<()> {
        // Sample paired-end FASTQ records (matching IDs)
        let read1 = b"@SEQ_ID1\nACGT\n+\n!!!!\n@SEQ_ID2\nTGCA\n+\n####\n";
        let read2 = b"@SEQ_ID1\nTTAA\n+\n$$$$\n@SEQ_ID2\nAATT\n+\n%%%%\n";

        // Create temp directory and files
        let tmp = tempdir()?;
        let out1_path = tmp.path().join("out1.fq.gz");
        let writer1 = fastq_writer(&out1_path, 10, 4, None).unwrap();
        let out2_path = tmp.path().join("out2.fq.gz");
        let writer2 = fastq_writer(&out2_path, 10, 4, None).unwrap();

        // No-op transformation
        let ranges: SortedSeqRanges = vec![SeqRange::From(3), SeqRange::To(2)]
            .into_iter()
            .collect();
        let mut actions = SubseqActions::builder();
        actions.add_action(SeqAction::Embed("UMI".to_string()), ranges);
        let actions = actions.build();
        let actions = SubseqPairedActions::new(Some(actions), None);
        // Run paired reader pipeline
        reader_seq_refine_paired_read(
            &mut read1.as_slice(),
            &mut Some(writer1),
            &mut read2.as_slice(),
            &mut Some(writer2),
            &actions,
            1,       // chunk_size
            10,      // buffer_size
            1,       // batch_size
            Some(1), // queue size
            1,
        )?;

        // Decompress and read output
        let mut buf1 = String::new();
        GzDecoder::new(File::open(out1_path)?).read_to_string(&mut buf1)?;

        let mut buf2 = String::new();
        GzDecoder::new(File::open(out2_path)?).read_to_string(&mut buf2)?;

        // Assert equality with original input
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
