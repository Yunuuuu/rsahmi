use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, Result};
use rustc_hash::FxHashSet as HashSet;

mod io;
mod mmap;
pub mod range;

use mmap::{
    mmap_kractor_paired_read, mmap_kractor_single_read,
    mmap_kractor_ubread_read,
};
use range::RangeKind;

pub fn mmap_kractor_reads(
    id_sets: HashSet<&[u8]>,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    ubread: Option<&str>,
    umi_pattern: Option<Vec<RangeKind>>,
    barcode_pattern: Option<Vec<RangeKind>>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    match (fq2, ubread) {
        (None, None) => mmap_kractor_single_read(
            id_sets,
            fq1,
            ofile1,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
        (Some(fq2), None) => {
            let ofile2 = ofile2.ok_or(anyhow!(
                "`ofile2` must be provided when processing paired-end reads"
            ))?;
            mmap_kractor_paired_read(
                id_sets,
                fq1,
                ofile1,
                fq2,
                ofile2,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
            )
        }
        (None, Some(ubread)) => {
            let umi_pattern = umi_pattern.ok_or(anyhow!(
                "`umi_pattern` must be provided when processing `ubread` reads"
            ))?;
            let barcode_pattern = barcode_pattern.ok_or(anyhow!(
                "`barcode_pattern` must be provided when processing `ubread` reads"
            ))?;
            mmap_kractor_ubread_read(
                id_sets,
                fq1,
                ofile1,
                ubread,
                umi_pattern,
                barcode_pattern,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
            )
        }
        (Some(_), Some(_)) => {
            return Err(anyhow!(
                "Both `fq2` and `ubread` cannot be provided simultaneously. Choose one."
            ));
        }
    }
}

// fn kractor_reads(
//     id_sets: HashSet<&[u8]>,
//     fq1: &str,
//     ofile1: &str,
//     fq2: Option<&str>,
//     ofile2: Option<&str>,
//     ubread: Option<&str>,
//     ub_pattern: Option<Vec<u8>>,
//     read_buffer: usize,
//     write_buffer: usize,
//     batch_size: usize,
//     read_queue: usize,
//     write_queue: usize,
//     threads: usize,
// ) -> std::result::Result<(), String> {
//     let (parser_tx, ref parser_rx) = bounded::<
//         Vec<(Vec<u8>, Option<Vec<u8>>, Option<Vec<u8>>)>,
//     >(threads * read_queue);
//     let (writer_tx, ref writer_rx) =
//         bounded::<Vec<Vec<u8>>>(threads * write_queue);

//     // reading lines
//     let mut fq1_reader = File::open(fq1).map_or_else(
//         |e| Err(format!("Cannot open {}: {}", fq1, e)),
//         |f| Ok(BufReader::with_capacity(read_buffer, f)),
//     )?;
//     let mut fq2_lines = fq2
//         .map(|file| {
//             File::open(file).map_err(|e| format!("Cannot open {}: {}", file, e))
//         })
//         .transpose()?
//         .map(|f| BufReader::with_capacity(read_buffer, f).split(b'\n'));
//     let mut ub_lines = ubread
//         .map(|file| {
//             File::open(file).map_err(|e| format!("Cannot open {}: {}", file, e))
//         })
//         .transpose()?
//         .map(|f| BufReader::with_capacity(read_buffer, f).split(b'\n'));
//     let buffer = fq1_reader
//         .fill_buf()
//         .map_err(|e| format!("Cannot read from fq1: {}", e))?;
//     let reader_batch_size = buffer.iter().filter(|b| **b == b'\n').count() + 1;
//     let mut fq1_lines = fq1_reader.split(b'\n');
//     let mut reader_batch = Vec::with_capacity(reader_batch_size);
//     while let Some(fq1_result) = fq1_lines.next() {
//         let fq1_line = fq1_result
//             .map_err(|e| format!("Failed to read bytes from {}: {}", fq1, e))?;
//         let fq2_line = if let Some(ref mut iter) = fq2_lines {
//             iter.next().transpose().map_err(|e| {
//                 format!(
//                     "Failed to read bytes from {}: {}",
//                     unsafe { fq2.unwrap_unchecked() },
//                     e
//                 )
//             })?
//         } else {
//             None
//         };
//         let ub_line = if let Some(ref mut iter) = ub_lines {
//             iter.next().transpose().map_err(|e| {
//                 format!(
//                     "Failed to read bytes from {}: {}",
//                     unsafe { ubread.unwrap_unchecked() },
//                     e
//                 )
//             })?
//         } else {
//             None
//         };
//         reader_batch.push((fq1_line, fq2_line, ub_line));
//         if reader_batch.len() >= reader_batch_size {
//             parser_tx
//                 .send(std::mem::take(&mut reader_batch))
//                 .map_err(|e| e.to_string())?;
//         }
//     }
//     // Send remaining records
//     if !reader_batch.is_empty() {
//         parser_tx.send(reader_batch).map_err(|e| e.to_string())?;
//     }
//     Ok(())
// }

pub fn read_sequence_id_from_koutput<P>(
    file: P,
    buffersize: usize,
) -> std::result::Result<Vec<Vec<u8>>, String>
where
    P: AsRef<Path> + Display,
{
    let opened =
        File::open(file).map_err(|e| format!("Open file failed: {}", e))?;
    let buffer = BufReader::with_capacity(buffersize, opened);
    let id_sets = buffer
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
        .collect::<Vec<Vec<u8>>>();
    Ok(id_sets)
}
