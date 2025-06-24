use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, Result};
use rustc_hash::FxHashSet as HashSet;

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
    umi_ranges: Option<Vec<RangeKind>>,
    barcode_ranges: Option<Vec<RangeKind>>,
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
            let umi_ranges = umi_ranges.ok_or(anyhow!(
                "`umi_ranges` must be provided when processing `ubread` reads"
            ))?;
            let barcode_ranges = barcode_ranges.ok_or(anyhow!(
                "`barcode_ranges` must be provided when processing `ubread` reads"
            ))?;
            mmap_kractor_ubread_read(
                id_sets,
                fq1,
                ofile1,
                ubread,
                umi_ranges,
                barcode_ranges,
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
