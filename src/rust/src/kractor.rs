use std::collections::HashSet;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use crossbeam_channel::bounded;
use extendr_api::prelude::*;
use memchr::memmem;
use noodles_fasta::io::Writer;
use noodles_fasta::record::{definition, sequence, Record};
use noodles_fastq::io::Reader;

#[extendr]
#[allow(clippy::too_many_arguments)]
fn kractor(
    koutput: &str,
    patterns: Robj,
    ofile: &str,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    io_buffer: usize,
    buffersize: usize,
    batchsize: usize,
    queue_capacity: usize,
    threads: usize,
) -> std::result::Result<(), String> {
    let patterns: Vec<&str> = patterns
        .as_str_vector()
        .ok_or("`patterns` must be a character vector")?;
    write_matching_output(
        koutput,
        &patterns,
        ofile,
        io_buffer,
        buffersize,
        batchsize,
        queue_capacity,
        threads,
    )?;
    let id_set = read_sequence_id_from_koutput(ofile, io_buffer)?;
    write_matching_reads(fq1, ofile1, fq2, ofile2, &id_set, io_buffer)
}

#[allow(clippy::too_many_arguments)]
fn write_matching_output<P>(
    koutput: P,
    patterns: &[&str],
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
    let patterns: Vec<_> = patterns
        .iter()
        .map(|s| memmem::Finder::new(s.as_bytes()))
        .collect();
    let needles = &patterns;

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
                            // SAFETY: pos is always in chunk
                            let line: Vec<u8> = chunk.drain(..= pos).collect();
                            let mut fields = line.splitn(3, |b| *b == b'\t');
                            if let Some(taxid) = fields.nth(2) {
                                if needles
                                    .iter()
                                    .any(|needle| needle.find(taxid).is_some())
                                {
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

fn read_sequence_id_from_koutput<P>(
    file: P,
    buffersize: usize,
) -> std::result::Result<HashSet<Vec<u8>>, String>
where
    P: AsRef<Path> + Display,
{
    // Collect IDs into a HashSet for fast lookup
    rprintln!("Extracting sequence IDs");
    let opened =
        File::open(file).map_err(|e| format!("Open file failed: {}", e))?;
    let buffer = BufReader::with_capacity(buffersize, opened);
    let id_set = buffer
        .lines()
        .filter_map(|line| {
            line.ok().and_then(|str| {
                // we selected the second column
                str.split("\t").nth(1).and_then(|second| {
                    // we remove empty sequence IDs
                    if second.is_empty() {
                        None
                    } else {
                        Some(second.as_bytes().to_vec())
                    }
                })
            })
        })
        .collect::<HashSet<Vec<u8>>>();
    Ok(id_set)
}

fn write_matching_reads<P>(
    fq1: P,
    ofile1: P,
    fq2: Option<P>,
    ofile2: Option<P>,
    id_set: &HashSet<Vec<u8>>,
    buffersize: usize,
) -> std::result::Result<(), String>
where
    P: AsRef<Path> + Display,
{
    // Iterate all FASTQ records
    rprintln!("Extracting the matching sequence from: {}", fq1);
    write_matching_records(fq1, ofile1, buffersize, id_set)?;

    if let (Some(in_file), Some(out_file)) = (fq2, ofile2) {
        rprintln!("Extracting the matching sequence from: {}", in_file);
        write_matching_records(in_file, out_file, buffersize, id_set)?;
    }
    Ok(())
}

fn write_matching_records<P>(
    fastq: P,
    fasta: P,
    buffersize: usize,
    id_set: &HashSet<Vec<u8>>,
) -> std::result::Result<(), String>
where
    P: AsRef<Path> + Display,
{
    // Open input FASTQ
    let in_file =
        File::open(fastq).map_err(|e| format!("Open file failed: {}", e))?;
    let in_buf = BufReader::with_capacity(buffersize, in_file);
    let mut reader = Reader::new(in_buf);

    // Open output FASTA
    let out_file =
        File::create(fasta).map_err(|e| format!("Open file failed: {}", e))?;
    let out_buf = BufWriter::with_capacity(buffersize, out_file);
    let mut writer = Writer::new(out_buf);

    // Iterate all FASTQ records
    for read in reader.records() {
        let record =
            read.map_err(|e| format!("Parse fastq record failed: {}", e))?;
        if id_set.contains(<&[u8]>::from(record.name())) {
            // convert FASTQ to FASTA record and write
            let definition = definition::Definition::new(
                record.name(),
                Some(record.description().to_owned()),
            );
            let sequence =
                sequence::Sequence::from(record.sequence().to_owned());
            let fasta_record = Record::new(definition, sequence);
            writer
                .write_record(&fasta_record)
                .map_err(|e| format!("Write fastq record failed: {}", e))?;
        }
    }
    writer
        .get_mut()
        .flush()
        .map_err(|e| format!("Flush lines failed: {}", e))?;
    Ok(())
}

extendr_module! {
    mod kractor;
    fn kractor;
}
