use std::io::BufReader;
use std::path::Path;

use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;

// use libdeflater::Compressor;
use crate::batchsender::BatchSender;
use crate::fastq_reader::*;
use crate::parser::fastq::FastqRecord;
use crate::seq_action::*;

pub(crate) fn reader_seq_refine_single_read<P: AsRef<Path> + ?Sized>(
    input_path: &P,
    input_bar: Option<ProgressBar>,
    output_path: &P,
    output_bar: Option<ProgressBar>,
    actions: &SubseqActions,
    compression_level: u32,
    chunk_size: usize,
    buffer_size: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let input: &Path = input_path.as_ref();
    let output: &Path = output_path.as_ref();
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
            let mut writer = fastq_writer(output, buffer_size, compression_level, output_bar)?;

            // Iterate over each received batch of records
            for chunk in writer_rx {
                for record in chunk {
                    record.write(&mut writer)?;
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
                let mut thread_tx = BatchSender::with_capacity(chunk_size, tx);
                while let Ok(records) = rx.recv() {
                    for mut record in records {
                        actions.transform_fastq_bytes(&mut record)?;
                        thread_tx.send(record).map_err(|e| {
                            anyhow!(
                                "(Parser) Failed to send parsed record to Writer thread: {}",
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
            let reader = fastq_reader(input, input_bar)?;
            let mut reader = FastqReader::new(BufReader::with_capacity(buffer_size, reader));
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
    use std::io::Write;

    use tempfile::NamedTempFile;

    use super::*;
    use crate::seq_action::{SeqAction, SeqRange, SortedSeqRanges, SubseqActions};

    fn dummy_fastq() -> &'static [u8] {
        b"@SEQ_ID\nGATTTGGGG\n+\nIIIIIIIII\n@SEQ_ID2\nTTACAGGGA\n+\nIIIIIIIII\n"
    }

    #[test]
    fn test_reader_seq_refine_single_read() {
        // Write dummy FASTQ to temp input file
        let mut input_file = NamedTempFile::new().expect("failed to create temp input file");
        input_file
            .write_all(dummy_fastq())
            .expect("failed to write FASTQ");

        let input_path = input_file.path().to_path_buf();

        // Prepare temp output file
        let output_file = NamedTempFile::new().expect("failed to create temp output file");
        let output_path = output_file.path().to_path_buf();

        // Define UMI extraction action
        let ranges: SortedSeqRanges = vec![SeqRange::From(3), SeqRange::To(2)]
            .into_iter()
            .collect();
        let mut builder = SubseqActions::builder();
        builder.add_action(SeqAction::Embed("UMI".to_string()), ranges);
        let actions = builder.build();

        // Call the function
        let result = reader_seq_refine_single_read(
            &input_path,
            None, // No progress bar
            &output_path,
            None, // No progress bar
            &actions,
            0,    // No compression
            1,       // chunk size
            8192,    // buffer size
            Some(2), // queue size
            2,       // threads
        );

        assert!(
            result.is_ok(),
            "reader_seq_refine_single_read failed: {:?}",
            result
        );

        // Check output contents
        let output_contents = std::fs::read_to_string(output_path).expect("failed to read output");
        assert!(
            output_contents.contains("@SEQ_ID RSAHMI{UMI:GATTGGGG}"),
            "Unexpected output: {}",
            output_contents
        );
        assert!(
            output_contents.contains("@SEQ_ID2 RSAHMI{UMI:TTCAGGGA}"),
            "Unexpected output: {}",
            output_contents
        );
    }
}
