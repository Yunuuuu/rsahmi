use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::iter::zip;

use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::batchsender::BatchSender;
use crate::parser::fastq::FastqRecord;
use crate::{fastq_reader, seq_action::*};

pub(crate) fn reader_seq_refine_paired_read<R1: Read + Send, R2: Read + Send>(
    reader1: R1,
    ofile1: Option<&str>,
    reader2: R2,
    ofile2: Option<&str>,
    ref actions: SubseqPairedActions,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Create output file and wrap in buffered writer
    let writer1;
    if let Some(file) = ofile1 {
        let file = File::create(&file)?;
        let bw = BufWriter::with_capacity(buffer_size, file);
        writer1 = Some(GzEncoder::new(bw, Compression::new(4)));
        // writer1 = Some(bw);
    } else {
        writer1 = None
    }
    let writer2;
    if let Some(file) = ofile2 {
        let file = File::create(&file)?;
        let bw = BufWriter::with_capacity(buffer_size, file);
        writer2 = Some(GzEncoder::new(bw, Compression::new(4)));
        // writer2 = Some(bw);
    } else {
        writer2 = None
    }

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<(FastqRecord<Bytes>, FastqRecord<Bytes>)>>,
            Receiver<Vec<(FastqRecord<Bytes>, FastqRecord<Bytes>)>>,
        ) = crate::new_channel(nqueue);

        let (writer1_tx, writer1_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let (writer2_tx, writer2_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let (reader1_tx, reader1_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let (reader2_tx, reader2_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let has_writer1 = writer1.is_some();
        let writer1_handle = if let Some(mut writer) = writer1 {
            Some(scope.spawn(move || -> Result<()> {
                for chunk in writer1_rx {
                    for record in chunk {
                        record.write(&mut writer)?;
                    }
                }
                writer.flush()?;
                Ok(())
            }))
        } else {
            None
        };

        let has_writer2 = writer2.is_some();
        let writer2_handle = if let Some(mut writer) = writer2 {
            Some(scope.spawn(move || -> Result<()> {
                for chunk in writer2_rx {
                    for record in chunk {
                        record.write(&mut writer)?;
                    }
                }
                writer.flush()?;
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
                    writer1_tx.send(records1)?;
                }
                if has_writer2 {
                    writer2_tx.send(records2)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let mut parser_handles = Vec::with_capacity(10);
        for _ in 0 .. 10 {
            let rx1 = reader1_rx.clone();
            let rx2 = reader2_rx.clone();
            let tx = parser_tx.clone();
            let handle =
                scope.spawn(move || -> Result<()> {
                    let mut thread_tx = BatchSender::with_capacity(batch_size, tx);
                    loop {
                        let (records1, records2) = match (rx1.recv(), rx2.recv()) {
                            (Ok(rec1), Ok(rec2)) => (rec1, rec2),
                            (Err(_), Ok(_)) => {
                                return Err(anyhow!(
                                    "FASTQ pairing error: read1 channel closed unexpectedly before read2"
                                ));
                            }
                            (Ok(_), Err(_)) => {
                                return Err(anyhow!(
                                    "FASTQ pairing error: read2 channel closed unexpectedly before read1"
                                ));
                            }
                            (Err(_), Err(_)) => {
                                return Ok(());
                            }
                        };
                        if records1.len() != records2.len() {
                            return Err(anyhow!("FASTQ pairing error: mismatched record counts"));
                        }
                        // Initialize a thread-local batch sender for matching records
                        for (mut record1, mut record2) in zip(records1, records2) {
                            if record1.id != record2.id {
                                return Err(anyhow!(
                                    "FASTQ pairing error: mismatched record id (read1: {}, read2: {})",
                                    String::from_utf8_lossy(&record1.id),
                                    String::from_utf8_lossy(&record2.id)
                                ));
                            }
                            actions.transform_fastq_bytes(&mut record1, &mut record2)?;
                            thread_tx.send((record1, record2))?;
                        }
                    }
                });
            parser_handles.push(handle);
        }
        drop(reader1_rx);
        drop(reader2_rx);
        drop(parser_tx);

        // ─── reader Thread ─────────────────────────────────────
        let reader1_handle = scope.spawn(move || -> Result<()> {
            let mut reader1 = fastq_reader::FastqReader::new(BufReader::new(reader1));
            let mut reader1_tx = BatchSender::with_capacity(chunk_size, reader1_tx);
            while let Some(record) = reader1.read_record()? {
                reader1_tx.send(record)?;
            }
            reader1_tx.flush()?;
            Ok(())
        });

        let reader2_handle = scope.spawn(move || -> Result<()> {
            let mut reader2 = fastq_reader::FastqReader::new(BufReader::new(reader2));
            let mut reader2_tx = BatchSender::with_capacity(chunk_size, reader2_tx);
            while let Some(record) = reader2.read_record()? {
                reader2_tx.send(record)?;
            }
            reader2_tx.flush()?;
            Ok(())
        });

        // ─── Join Threads and Propagate Errors ────────────────
        if let Some(writer_handle) = writer1_handle {
            writer_handle
                .join()
                .map_err(|e| anyhow!("Writer1 thread panicked: {:?}", e))??
        };
        if let Some(writer_handle) = writer2_handle {
            writer_handle
                .join()
                .map_err(|e| anyhow!("Writer2 thread panicked: {:?}", e))??
        };
        writer_handle
            .join()
            .map_err(|e| anyhow!("Writer thread panicked: {:?}", e))??;

        // we ensure no data to send
        for handler in parser_handles {
            handler
                .join()
                .map_err(|e| anyhow!("Parser thread panicked: {:?}", e))??;
        }
        reader1_handle
            .join()
            .map_err(|e| anyhow!("Reader1 thread panicked: {:?}", e))??;
        reader2_handle
            .join()
            .map_err(|e| anyhow!("Reader2 thread panicked: {:?}", e))??;
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

    #[test]
    fn test_reader_seq_refine_paired_read_passthrough() -> Result<()> {
        // Sample paired-end FASTQ records (matching IDs)
        let read1 = b"@SEQ_ID1\nACGT\n+\n!!!!\n@SEQ_ID2\nTGCA\n+\n####\n";
        let read2 = b"@SEQ_ID1\nTTAA\n+\n$$$$\n@SEQ_ID2\nAATT\n+\n%%%%\n";

        // Create temp directory and files
        let tmp = tempdir()?;
        let out1_path = tmp.path().join("out1.fq.gz");
        let out2_path = tmp.path().join("out2.fq.gz");

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
            &read1[..],
            Some(out1_path.to_str().unwrap()),
            &read2[..],
            Some(out2_path.to_str().unwrap()),
            actions,
            1,       // chunk_size
            4096,    // buffer_size
            1,       // batch_size
            Some(4), // queue size
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
            b"@SEQ_ID1 RSAHMI{UMI:TTA}\nTTAA\n+\n$$$$\n@SEQ_ID2 RSAHMI{UMI:AAT}\nAATT\n+\n%%%%\n"
        );
        Ok(())
    }
}
