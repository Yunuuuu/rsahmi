use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use aho_corasick::AhoCorasick;
use crossbeam_channel::bounded;
use extendr_api::prelude::*;

#[allow(clippy::too_many_arguments)]
pub fn write_matching_output<P>(
    koutput: P,
    needles: &AhoCorasick,
    ofile: P,
    io_buffer: usize,
    buffersize: usize,
    batchsize: usize,
    queue_capacity: usize,
    threads: usize,
) -> std::result::Result<(), String>
where
    P: AsRef<Path> + Display,
{
    rprintln!("Extracting matching kraken2 output from: {}", koutput);
    // one thread is kept for writer
    let threads = if threads <= 2 { 1 } else { threads - 1 };
    let mut input = File::open(koutput).map_err(|e| e.to_string())?;
    let mut output = BufWriter::with_capacity(
        io_buffer,
        File::create(ofile).map_err(|e| e.to_string())?,
    );

    // Start the processor threads
    let (writer_tx, ref writer_rx) =
        bounded::<Vec<Vec<u8>>>(threads * queue_capacity);
    let (work_tx, ref work_rx) = bounded::<Vec<u8>>(threads * queue_capacity);

    std::thread::scope(|scope| {
        // let patterns = Arc::new(patterns);

        // Writer thread: write data to the file
        // accept lines from `writer_rx`
        let writer_handle =
            scope.spawn(move || -> std::result::Result<(), String> {
                for batch in writer_rx.iter() {
                    for line in batch {
                        output
                            .write_all(&line)
                            .map_err(|e| format!("Write failed: {e}"))?;
                    }
                }
                output.flush().map_err(|e| format!("Flush failed: {e}"))?;
                Ok(())
            });

        // Worker threads, will send each passed lines to `writer`
        let mut worker_handles = Vec::with_capacity(threads);
        for _ in 0 .. threads {
            let tx = writer_tx.clone();
            let handle =
                scope.spawn(move || -> std::result::Result<(), String> {
                    let mut batch_buffer = Vec::with_capacity(batchsize);
                    for mut chunk in work_rx.iter() {
                        while let Some(pos) =
                            chunk.iter().position(|&b| b == b'\n')
                        {
                            let line: Vec<u8> = chunk.drain(..= pos).collect();
                            let mut fields = line.splitn(3, |b| *b == b'\t');
                            if let Some(taxid) = fields.nth(2) {
                                if needles.find(taxid).is_some() {
                                    batch_buffer.push(line);
                                    if batch_buffer.len()
                                        == batch_buffer.capacity()
                                    {
                                        // Send batch by moving it, avoid cloning
                                        let send_batch =
                                            std::mem::take(&mut batch_buffer);
                                        tx.send(send_batch).map_err(|e| {
                                            format!(
                                                "Send to writer failed: {e}"
                                            )
                                        })?;
                                        batch_buffer =
                                            Vec::with_capacity(batchsize)
                                    }
                                }
                            }
                        }
                    }
                    if !batch_buffer.is_empty() {
                        tx.send(batch_buffer).map_err(|e| {
                            format!("Send to writer failed: {e}")
                        })?;
                    }
                    Ok(())
                });
            worker_handles.push(handle);
        }

        // read files and pass lines
        let mut buffer = vec![0; buffersize]; // 128 KB
        let mut offset: usize = 0;
        loop {
            let read_bytes = input
                .read(&mut buffer[offset ..])
                .map_err(|e| format!("Read failed: {e}"))?;
            if read_bytes == 0 {
                // ensure all buffer has been send
                if offset > 0 {
                    let mut chunk =
                        buffer.drain(..= offset).collect::<Vec<u8>>();
                    // always ensure the last line has `\n`
                    if !chunk.ends_with(b"\n") {
                        chunk.push(b'\n');
                    }
                    work_tx
                        .send(chunk)
                        .map_err(|e| format!("Send to workers failed: {e}"))?;
                }
                break;
            }
            offset += read_bytes;
            match buffer[.. offset].iter().rposition(|b| *b == b'\n') {
                Some(pos) => {
                    offset -= pos + 1;
                    let chunk = buffer.drain(..= pos).collect::<Vec<u8>>();
                    work_tx
                        .send(chunk)
                        .map_err(|e| format!("Send to workers failed: {e}"))?;
                    buffer.resize(buffer.capacity(), 0);
                }
                None => {
                    if buffer.capacity() < offset {
                        buffer
                            .try_reserve(
                                buffer.capacity() - buffer.len() + buffersize,
                            )
                            .map_err(|e| {
                                format!("Increase buffer failed: {e}")
                            })?;
                        buffer.resize(buffer.capacity(), 0);
                    }
                }
            };
        }

        // close work channel to stop workers
        drop(work_tx);

        // we ensure no lines to send
        for handler in worker_handles {
            let _ = handler.join().map_err(|_| "worker thread panicked")?;
        }
        drop(writer_tx); // close writer channel to stop writer

        // we ensure all lines have been writen
        let _ = writer_handle.join().map_err(|_| "writer thread panicked")?;
        Ok(())
    })
}
