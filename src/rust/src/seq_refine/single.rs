use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use libdeflater::{CompressionLvl, Compressor};

use super::seq_action::*;
use crate::batchsender::BatchSender;
use crate::fastq_reader::*;
use crate::fastq_record::FastqRecord;
use crate::utils::*;

pub(crate) fn seq_refine_single_read<P: AsRef<Path> + ?Sized>(
    input_path: &P,
    input_bar: Option<ProgressBar>,
    output_path: &P,
    output_bar: Option<ProgressBar>,
    actions: &SubseqActions,
    compression_level: i32,
    batch_size: usize,
    chunk_bytes: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let input: &Path = input_path.as_ref();
    let output: &Path = output_path.as_ref();

    // Ensure compression level is validated and converted before entering thread scope.
    // Doing this outside avoids redundant validation across parser threads.
    let compression_level = CompressionLvl::new(compression_level)
        .map_err(|e| anyhow!("Invalid 'compression_level': {:?}", e))?;
    std::thread::scope(|scope| -> Result<()> {
        // Two communication pipelines are set up to decouple IO and CPU-intensive work:
        // - reader_tx: transfers raw FASTQ records to parser threads
        // - writer_tx: receives compressed byte chunks from parser threads
        let (writer_tx, writer_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = new_channel(nqueue);
        let (reader_tx, reader_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // A single thread handles file output to ensure atomic write order and leverage buffered IO.
        // This thread consumes compressed chunks, not raw records, for performance.
        let writer_handle = scope.spawn(move || -> Result<()> {
            let mut writer = BufWriter::with_capacity(chunk_bytes, new_writer(output, output_bar)?);

            // Iterate over each received batch of records
            for chunk in writer_rx {
                writer
                    .write_all(&chunk)
                    .with_context(|| format!("(Writer) Failed to write FastqRecord to output"))?;
            }
            writer
                .flush()
                .with_context(|| format!("(Writer) Failed to flush writer"))?;
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Multiple parser threads are used to exploit CPU parallelism for `transform_fastq` and gzip compression.
        // Each thread transforms records and buffers them into a local pool,
        // which is periodically flushed into the writer pipeline.
        let mut parser_handles = Vec::with_capacity(threads);
        let gzip = gz_compressed(output);
        for _ in 0 .. threads {
            let rx = reader_rx.clone();
            let tx = writer_tx.clone();
            let handle = scope.spawn(move || -> Result<()> {
                // Temporary buffer for current output chunk
                let mut records_pool: Vec<u8> = Vec::with_capacity(chunk_bytes);
                let mut compressor = Compressor::new(compression_level);
                while let Ok(records) = rx.recv() {
                    for mut record in records {
                        // Apply trimming, tag embedding, and other sequence transformations
                        actions.transform_fastq(&mut record)?;

                        // Flush when pool is too full to accept the next record.
                        // This ensures output chunks remain near the target block size.
                        if records_pool.capacity() - records_pool.len() < record.bytes_size() {
                            let mut pack = Vec::with_capacity(chunk_bytes);
                            std::mem::swap(&mut records_pool, &mut pack);
                            // Compress if gzip file
                            if gzip {
                                pack = gzip_pack(&pack, &mut compressor)?
                            }

                            // Send compressed or raw bytes to writer
                            tx.send(pack).with_context(|| {
                                format!("(Parser) Failed to send parsed record to Writer thread")
                            })?;
                        }

                        // Append encoded record to buffer
                        record.extend(&mut records_pool);
                    }
                }

                // Flush remaining records if any
                if !records_pool.is_empty() {
                    let pack = if gzip {
                        gzip_pack(&records_pool, &mut compressor)?
                    } else {
                        records_pool
                    };
                    tx.send(pack).with_context(|| {
                        format!("(Parser) Failed to send parsed record to Writer thread")
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
            let mut reader =
                FastqReader::with_capacity(BUFFER_SIZE, new_reader(input, BUFFER_SIZE, input_bar)?);
            let mut reader_tx = BatchSender::with_capacity(batch_size, reader_tx);
            while let Some(record) = reader
                .read_record()
                .with_context(|| format!("(Reader) Failed to read FASTQ record"))?
            {
                reader_tx.send(record).with_context(|| {
                    format!("(Reader) Failed to send FASTQ records to Parser thread")
                })?;
            }
            reader_tx.flush().with_context(|| {
                format!("(Reader) Failed to flush FASTQ records to Parser thread")
            })?;
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
    use super::{SeqAction, SubseqActions};
    use crate::seq_range::{SeqRange, SeqRanges};

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
        let ranges: SeqRanges = vec![SeqRange::From(3), SeqRange::To(2)]
            .into_iter()
            .collect();
        let mut builder = SubseqActions::builder();
        builder
            .add_action(
                SeqAction::Embed(Bytes::from_owner("UMI".as_bytes())),
                ranges,
            )
            .unwrap();
        let actions = builder.build().unwrap();

        // Call the function
        let result = seq_refine_single_read(
            &input_path,
            None, // No progress bar
            &output_path,
            None, // No progress bar
            &actions,
            1,       // No compression
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
