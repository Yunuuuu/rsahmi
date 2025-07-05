use std::fs::File;

use anyhow::{anyhow, Result};
use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};
use memmap2::Mmap;

mod paired;
mod single;

use crate::reader::slice::SliceProgressBarReader;
use crate::seq_action::*;

pub(crate) fn mmap_seq_refine(
    fq1: &str,
    ofile1: Option<&str>,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    actions1: Option<SubseqActions>,
    actions2: Option<SubseqActions>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    if let Some(fq2) = fq2 {
        mmap_seq_refine_paired_read(
            fq1,
            ofile1,
            fq2,
            ofile2,
            actions1,
            actions2,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        )
    } else {
        mmap_seq_refine_single_read(
            fq1,
            ofile1,
            actions1,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        )
    }
}

fn mmap_seq_refine_single_read(
    fq1: &str,
    ofile1: Option<&str>,
    actions: Option<SubseqActions>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    let ofile1 = ofile1.ok_or_else(|| anyhow!("No output file specified."))?;
    let actions = actions.ok_or_else(|| anyhow!("No sequence actions were specified."))?;
    let file = File::open(fq1)?;
    let map = unsafe { Mmap::map(&file) }?;
    crate::mmap_advice(&map)?;

    let mut reader = SliceProgressBarReader::new(&map);
    reader.set_label("fq1");
    let progress_style = ProgressStyle::with_template(
        "{prefix:.bold.green} {wide_bar:.cyan/blue} {decimal_bytes}/{decimal_total_bytes} [{elapsed_precise}] {decimal_bytes_per_sec} ({eta})",
    )?;
    let pb = ProgressBar::new(map.len() as u64).with_finish(ProgressFinish::Abandon);
    pb.set_prefix("Parsing fq1");
    pb.set_style(progress_style);

    #[cfg(not(test))]
    reader.attach_bar(pb);

    single::mmap_seq_refine_single_read(
        reader,
        ofile1,
        actions,
        chunk_size,
        buffer_size,
        batch_size,
        nqueue,
    )
}

fn mmap_seq_refine_paired_read(
    fq1: &str,
    ofile1: Option<&str>,
    fq2: &str,
    ofile2: Option<&str>,
    actions1: Option<SubseqActions>,
    actions2: Option<SubseqActions>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    if ofile1.is_none() && ofile2.is_none() {
        return Err(anyhow!("No output file specified."));
    }
    if actions1.is_none() && actions2.is_none() {
        return Err(anyhow!(
            "No sequence actions were specified. Please provide at least one action to proceed"
        ));
    }
    let file1 = File::open(fq1)?;
    let map1 = unsafe { Mmap::map(&file1) }?;
    crate::mmap_advice(&map1)?;

    let file2 = File::open(fq2)?;
    let map2 = unsafe { Mmap::map(&file2) }?;
    crate::mmap_advice(&map2)?;

    let mut reader1 = SliceProgressBarReader::new(&map1);
    reader1.set_label("fq1");

    let mut reader2 = SliceProgressBarReader::new(&map2);
    reader2.set_label("fq2");

    let progress_style = ProgressStyle::with_template(
        "{prefix:.bold.green} {wide_bar:.cyan/blue} {decimal_bytes}/{decimal_total_bytes} [{elapsed_precise}] {decimal_bytes_per_sec} ({eta})",
    )?;
    let progress = MultiProgress::new();
    let pb1 =
        progress.add(ProgressBar::new(map1.len() as u64).with_finish(ProgressFinish::Abandon));
    pb1.set_prefix("Parsing fq1");
    pb1.set_style(progress_style.clone());

    #[cfg(not(test))]
    reader1.attach_bar(pb1);

    let pb2 =
        progress.add(ProgressBar::new(map2.len() as u64).with_finish(ProgressFinish::Abandon));
    pb2.set_prefix("Parsing fq2");
    pb2.set_style(progress_style);

    #[cfg(not(test))]
    reader2.attach_bar(pb2);

    let actions = SubseqPairedActions::new(actions1, actions2);
    paired::mmap_seq_refine_paired_read(
        reader1,
        ofile1,
        reader2,
        ofile2,
        actions,
        chunk_size,
        buffer_size,
        batch_size,
        nqueue,
    )
}
