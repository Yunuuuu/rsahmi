use std::fs::File;
use std::io::{BufWriter, Read};
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};
use std::sync::Arc;

use anyhow::{anyhow, Error, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};

use crate::batchsender::BatchSender;
use crate::parser::fastq::FastqBytesChunkReader;
use crate::parser::fastq::FastqContainer;
use crate::parser::fastq::FastqRecord;
use crate::reader::bytes::BytesReader;
use crate::seq_action::*;

pub(crate) fn reader_seq_refine_single_read<R: Read + Send>(
    reader: R,
    ofile: &str,
    ref actions: SubseqActions,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    // Open output file and wrap in buffered writer
    let mut writer = BufWriter::with_capacity(buffer_size, File::create(ofile)?);

    std::thread::scope(|scope| -> Result<()> {
        // Create a channel between the parser and writer threads
        // The channel transmits batches (Vec<FastqRecord>)
        let (parser_tx, writer_rx): (
            Sender<Vec<FastqRecord<Bytes>>>,
            Receiver<Vec<FastqRecord<Bytes>>>,
        ) = crate::new_channel(nqueue);

        // ─── Writer Thread ─────────────────────────────────────
        // Consumes batches of records and writes them to file
        let writer_handle = scope.spawn(move || -> Result<()> {
            // Iterate over each received batch of records
            for chunk in writer_rx {
                for record in chunk {
                    record.write(&mut writer)?;
                }
            }
            Ok(())
        });

        // ─── Parser Thread ─────────────────────────────────────
        // Streams FASTQ data, filters by ID set, sends batches to writer
        let mut reader = BytesReader::new(reader);
        reader.set_label("fq1");
        let mut reader = FastqBytesChunkReader::with_capacity(chunk_size, reader);

        let parser_handle = scope.spawn(move || -> Result<()> {
            // will move `reader`, `parser_tx`, and `id_sets`
            // Shared atomic flag to signal if any thread encountered an error
            let has_error = Arc::new(AtomicBool::new(false));
            // A bounded channel to capture the first error that occurs (capacity = 1)
            let (err_tx, err_rx) = crossbeam_channel::bounded(1);
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
                    let thread_err_tx: Sender<Error> = err_tx.clone();
                    s.spawn(move |_| {
                        let mut container = FastqContainer::new();
                        loop {
                            match chunk.read_record(&mut container) {
                                Ok(opt_record) => match opt_record {
                                    Some(mut record) => {
                                        match actions.transform_fastq_bytes(&mut record) {
                                            Ok(_) => {
                                                // Attempt to send the matching record to the writer thread.
                                                match thread_tx.send(record) {
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
                                            Err(e) => {
                                                // Mark that an error has occurred
                                                thread_has_error.store(true, Relaxed);
                                                // Only the first error is reported
                                                let _ = thread_err_tx.try_send(e.into());
                                                return ();
                                            }
                                        }
                                    }
                                    None => break,
                                },
                                Err(e) => {
                                    // Mark that an error has occurred
                                    thread_has_error.store(true, Relaxed);
                                    // Only the first error is reported
                                    let _ = thread_err_tx.try_send(e.into());
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
                Err(err)
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
    use std::fs;
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_reader_seq_refine_single_read_embed_trim() -> Result<()> {
        let fastq_data = b"@SEQ_ID\nACGTGATCGT\n+\nFFFFFFFFFF\n";
        let reader = Cursor::new(&fastq_data[..]);

        let tmpfile = tempfile::NamedTempFile::new()?;
        let out_path = tmpfile.path().to_str().unwrap();

        // Define a seq_range and embed_trim action
        let ranges: SortedSeqRanges = vec![SeqRange::From(8), SeqRange::To(2)]
            .into_iter()
            .collect();
        let mut actions = SubseqActions::builder();
        actions.add_action(SeqAction::Embed("UMI".to_string()), ranges);
        let actions = actions.build();

        reader_seq_refine_single_read(reader, out_path, actions, 1, 1024, 10, Some(2))?;

        let output = fs::read_to_string(out_path)?;
        assert!(output.contains("RSAHMI{UMI:ACGT}"));

        Ok(())
    }
}
