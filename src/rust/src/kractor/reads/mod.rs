use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, Result};
use rustc_hash::FxHashSet as HashSet;

mod paired;
mod single;

use indicatif::{MultiProgress, ProgressBar, ProgressFinish};

use crate::utils::*;

pub(super) fn kractor_reads(
    koutput: &str,
    fq1: &str,
    ofile1: Option<&str>,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    compression_level: i32,
    batch_size: usize,
    chunk_bytes: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let ids = read_sequence_id_from_koutput(koutput, 126 * 1024)
        .map_err(|e| anyhow!("Failed to read sequence IDs: {}", e))?;
    let id_sets = ids
        .iter()
        .map(|id| id.as_slice())
        .collect::<HashSet<&[u8]>>();
    let threads = threads.max(1); // always use at least one thread
    if let Some(fq2) = fq2 {
        kractor_reads_paired(
            &id_sets,
            fq1,
            ofile1,
            fq2,
            ofile2,
            batch_size,
            chunk_bytes,
            compression_level,
            nqueue,
            threads,
        )
    } else {
        kractor_reads_single(
            &id_sets,
            fq1,
            ofile1,
            batch_size,
            chunk_bytes,
            compression_level,
            nqueue,
            threads,
        )
    }
}

fn kractor_reads_single(
    id_sets: &HashSet<&[u8]>,
    fq1: &str,
    ofile1: Option<&str>,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: i32,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let ofile1 = ofile1.ok_or_else(|| anyhow!("No output file specified."))?;
    let reader_style = progress_reader_style()?;
    let writer_style = progress_writer_style()?;
    let progress = MultiProgress::new();
    let pb1 = progress.add(
        ProgressBar::new(std::fs::metadata(fq1)?.len() as u64).with_finish(ProgressFinish::Abandon),
    );
    pb1.set_prefix("Reading fastq");
    pb1.set_style(reader_style);

    let pb2 = progress.add(ProgressBar::no_length().with_finish(ProgressFinish::Abandon));
    pb2.set_prefix("Writing fastq");
    pb2.set_style(writer_style);

    single::parse_single(
        id_sets,
        &fq1,
        Some(pb1),
        &ofile1,
        Some(pb2),
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    )
}

fn kractor_reads_paired(
    id_sets: &HashSet<&[u8]>,
    fq1: &str,
    ofile1: Option<&str>,
    fq2: &str,
    ofile2: Option<&str>,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: i32,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    if ofile1.is_none() && ofile2.is_none() {
        return Err(anyhow!("No output file specified."));
    }

    let reader_style = progress_reader_style()?;
    let writer_style = progress_writer_style()?;
    let progress = MultiProgress::new();
    let pb1 = progress.add(
        ProgressBar::new(std::fs::metadata(fq1)?.len() as u64).with_finish(ProgressFinish::Abandon),
    );
    pb1.set_prefix("Reading fq1");
    pb1.set_style(reader_style.clone());
    let pb2 = if let Some(_) = ofile1 {
        let pb2 = progress.add(ProgressBar::no_length().with_finish(ProgressFinish::Abandon));
        pb2.set_prefix("Writing fq1");
        pb2.set_style(writer_style.clone());
        Some(pb2)
    } else {
        None
    };

    let pb3 = progress.add(
        ProgressBar::new(std::fs::metadata(fq2)?.len() as u64).with_finish(ProgressFinish::Abandon),
    );
    pb3.set_prefix("Reading fq2");
    pb3.set_style(reader_style);
    let pb4 = if let Some(_) = ofile2 {
        let pb4 = progress.add(ProgressBar::no_length().with_finish(ProgressFinish::Abandon));
        pb4.set_prefix("Writing fq2");
        pb4.set_style(writer_style);
        Some(pb4)
    } else {
        None
    };
    paired::parse_paired(
        id_sets,
        fq1,
        Some(pb1),
        fq2,
        Some(pb3),
        ofile1,
        pb2,
        ofile2,
        pb4,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    )
}

fn read_sequence_id_from_koutput<P>(
    file: P,
    buffersize: usize,
) -> std::result::Result<Vec<Vec<u8>>, String>
where
    P: AsRef<Path> + Display,
{
    let opened = File::open(file).map_err(|e| format!("Open file failed: {}", e))?;
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
