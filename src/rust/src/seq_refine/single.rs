use std::io::{BufReader, Read, Write};

use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};

use crate::batchsender::BatchSender;
use crate::fastq_reader::*;
use crate::parser::fastq::FastqRecord;
use crate::seq_action::*;

pub(crate) fn reader_seq_refine_single_read<R: Read + Send, W: Write + Send>(
    reader: &mut R,
    writer: &mut W,
    actions: &SubseqActions,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let mut reader = FastqReader::new(BufReader::with_capacity(buffer_size, reader));
    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (writer_tx, writer_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        let (reader_tx, reader_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for record in chunk {
                    record.write(writer)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        let mut parser_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                let mut thread_tx = BatchSender::with_capacity(batch_size, tx);
                while let Ok(records) = rx.recv() {
                    for mut record in records {
                        actions.transform_fastq_bytes(&mut record)?;
                        thread_tx.send(record).map_err(|e| {
                            anyhow!(
                                "(Parser) Failed to send parsed record pair to Writer thread: {}",
                                e
                            )
                        })?;
                    }
                }
                thread_tx.flush().map_err(|e| {
                    anyhow!(
                        "(Parser) Failed to flush remaining records to Writer thread: {}",
                        e
                    )
                })?;
                Ok(())
            });
            parser_handles.push(handle);
        }
        drop(reader_rx);
        drop(writer_tx);

        // ─── reader Thread ─────────────────────────────────────
        let reader_handle = scope.spawn(move || -> Result<()> {
            let mut reader_tx = BatchSender::with_capacity(chunk_size, reader_tx);
            while let Some(record) = reader
                .read_record()
                .map_err(|e| anyhow!("(Reader) Error while reading FASTQ record: {}", e))?
            {
                reader_tx.send(record).map_err(|e| {
                    anyhow!(
                        "(Reader) Failed to send FASTQ record to Parser thread: {}",
                        e
                    )
                })?;
            }
            reader_tx
                .flush()
                .map_err(|e| anyhow!("(Reader) Failed to flush records to Parser thread: {}", e))?;
            Ok(())
        });

        // ─── Join Threads and Propagate Errors ────────────────
        writer_handle
            .join()
            .map_err(|e| anyhow!("(Writer) thread panicked: {:?}", e))??;
        for handler in parser_handles {
            handler
                .join()
                .map_err(|e| anyhow!("(Parser) thread panicked: {:?}", e))??;
        }
        reader_handle
            .join()
            .map_err(|e| anyhow!("(Reader) thread panicked: {:?}", e))??;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn dummy_fastq() -> &'static [u8] {
        b"@SEQ_ID\nGATTTGGGG\n+\nIIIIIIIII\n@SEQ_ID2\nTTACAGGGA\n+\nIIIIIIIII\n"
    }

    #[test]
    fn test_reader_seq_refine_single_read() {
        let mut input = Cursor::new(dummy_fastq());

        let mut output = Vec::new();
        let ranges: SortedSeqRanges = vec![SeqRange::From(3), SeqRange::To(2)]
            .into_iter()
            .collect();
        let mut actions = SubseqActions::builder();
        actions.add_action(SeqAction::Embed("UMI".to_string()), ranges);
        let actions = actions.build();

        let result =
            reader_seq_refine_single_read(&mut input, &mut output, &actions, 1, 1, 10, Some(2), 2);

        assert!(
            result.is_ok(),
            "reader_seq_refine_single_read failed: {:?}",
            result
        );

        // Optionally, validate output content:
        let out_str = String::from_utf8_lossy(&output);
        assert!(out_str.contains("@SEQ_ID RSAHMI{UMI:GATTGGGG}"));
        assert!(out_str.contains("@SEQ_ID2 RSAHMI{UMI:TTCAGGGA}"));
    }
}
