use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{anyhow, Result};
use rustc_hash::FxHashSet as HashSet;

pub(crate) mod io;
pub(crate) mod mmap;

use indicatif::{MultiProgress, ProgressBar, ProgressFinish};
use io::{reader_kractor_paired_read, reader_kractor_single_read};
use memmap2::Mmap;
use mmap::{mmap_kractor_paired_read, mmap_kractor_single_read};

use crate::reader::bytes::BytesProgressBarReader;
use crate::reader::slice::SliceProgressBarReader;

pub(crate) fn kractor_reads(
    koutput: &str,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    mmap: bool,
) -> Result<()> {
    println!("Extracting sequence IDs");
    let ids = read_sequence_id_from_koutput(koutput, 126 * 1024)
        .map_err(|e| anyhow!("Failed to read sequence IDs: {}", e))?;
    let id_sets = ids
        .iter()
        .map(|id| id.as_slice())
        .collect::<HashSet<&[u8]>>();
    println!(
        "Extracting the matching sequence from: {} {}",
        fq1,
        fq2.map_or_else(|| String::from(""), |s| format!("and {}", s))
    );
    let file = File::open(fq1)?;
    let style = crate::progress_style()?;

    if let Some(fq2) = fq2 {
        let ofile2 = ofile2.ok_or(anyhow!(
            "`ofile2` must be provided when processing paired-end reads"
        ))?;
        let file2 = File::open(fq2)?;
        // add progres bar
        let progress = MultiProgress::new();
        if mmap {
            let map1 = unsafe { Mmap::map(&file) }?;
            crate::mmap_advice(&map1)?;
            let mut reader1 = SliceProgressBarReader::new(&map1);
            reader1.set_label("read1");

            let map2 = unsafe { Mmap::map(&file2) }?;
            crate::mmap_advice(&map2)?;
            let mut reader2 = SliceProgressBarReader::new(&map2);
            reader2.set_label("read2");

            let pb1 = progress
                .add(ProgressBar::new(map1.len() as u64).with_finish(ProgressFinish::Abandon));
            pb1.set_prefix("Parsing read1");
            pb1.set_style(style.clone());

            #[cfg(not(test))]
            reader1.attach_bar(pb1);

            let pb2 = progress
                .add(ProgressBar::new(map2.len() as u64).with_finish(ProgressFinish::Abandon));
            pb2.set_prefix("Parsing read2");
            pb2.set_style(style);

            #[cfg(not(test))]
            reader2.attach_bar(pb2);

            mmap_kractor_paired_read(
                id_sets,
                reader1,
                ofile1,
                reader2,
                ofile2,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
            )
        } else {
            let mut reader1 = BytesProgressBarReader::new(&file);
            reader1.set_label("read1");

            let mut reader2 = BytesProgressBarReader::new(&file2);
            reader2.set_label("read2");

            let pb1 = progress.add(
                ProgressBar::new(file.metadata()?.len() as u64)
                    .with_finish(ProgressFinish::Abandon),
            );
            pb1.set_prefix("Parsing read1");
            pb1.set_style(style.clone());

            #[cfg(not(test))]
            reader1.attach_bar(pb1);

            let pb2 = progress.add(
                ProgressBar::new(file2.metadata()?.len() as u64)
                    .with_finish(ProgressFinish::Abandon),
            );
            pb2.set_prefix("Parsing read2");
            pb2.set_style(style);

            #[cfg(not(test))]
            reader2.attach_bar(pb2);

            reader_kractor_paired_read(
                id_sets,
                reader1,
                ofile1,
                reader2,
                ofile2,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
            )
        }
    } else {
        if mmap {
            let map = unsafe { Mmap::map(&file) }?;
            crate::mmap_advice(&map)?;
            let mut reader = SliceProgressBarReader::new(&map);
            reader.set_label("reads");

            let pb = ProgressBar::new(map.len() as u64).with_finish(ProgressFinish::Abandon);
            pb.set_prefix("Parsing reads");
            pb.set_style(style);

            #[cfg(not(test))]
            reader.attach_bar(pb);

            mmap_kractor_single_read(
                id_sets,
                reader,
                ofile1,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
            )
        } else {
            let mut reader = BytesProgressBarReader::new(&file);
            reader.set_label("reads");

            let pb = ProgressBar::new(file.metadata()?.len() as u64)
                .with_finish(ProgressFinish::Abandon);
            pb.set_prefix("Parsing reads");
            pb.set_style(style);

            #[cfg(not(test))]
            reader.attach_bar(pb);

            reader_kractor_single_read(
                id_sets,
                reader,
                ofile1,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
            )
        }
    }
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
