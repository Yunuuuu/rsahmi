use std::fs::File;
use std::path::Path;

use anyhow::{anyhow, Result};
use flate2::read::MultiGzDecoder;
use indicatif::{MultiProgress, ProgressBar, ProgressFinish};

pub(crate) mod paired;
pub(crate) mod single;

use crate::reader::bytes::ProgressBarReader;
use crate::seq_action::*;

pub(crate) fn reader_seq_refine(
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
        reader_seq_refine_paired_read(
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
        reader_seq_refine_single_read(
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

fn reader_seq_refine_single_read(
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
    let path: &Path = fq1.as_ref();
    let file: File = File::open(path)?;

    let style = crate::progress_style()?;
    let pb = ProgressBar::new(file.metadata()?.len() as u64).with_finish(ProgressFinish::Abandon);
    pb.set_prefix("Parsing fq1");
    pb.set_style(style);
    let reader = ProgressBarReader::new(&file, pb);

    if crate::gz_compressed(path) {
        single::reader_seq_refine_single_read(
            MultiGzDecoder::new(reader),
            ofile1,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        )
    } else {
        single::reader_seq_refine_single_read(
            reader,
            ofile1,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        )
    }
}

fn reader_seq_refine_paired_read(
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
    let path1: &Path = fq1.as_ref();
    let file1 = File::open(path1)?;
    let path2: &Path = fq2.as_ref();
    let file2 = File::open(path2)?;

    let style = crate::progress_style()?;
    let progress = MultiProgress::new();
    let pb1 = progress
        .add(ProgressBar::new(file1.metadata()?.len() as u64).with_finish(ProgressFinish::Abandon));
    pb1.set_prefix("Parsing fq1");
    pb1.set_style(style.clone());
    let reader1 = ProgressBarReader::new(file1, pb1);

    let pb2 = progress
        .add(ProgressBar::new(file2.metadata()?.len() as u64).with_finish(ProgressFinish::Abandon));
    pb2.set_prefix("Parsing fq2");
    pb2.set_style(style);
    let reader2 = ProgressBarReader::new(file2, pb2);

    let actions = SubseqPairedActions::new(actions1, actions2);
    match (crate::gz_compressed(path1), crate::gz_compressed(path2)) {
        (true, true) => paired::reader_seq_refine_paired_read(
            MultiGzDecoder::new(reader1),
            ofile1,
            MultiGzDecoder::new(reader2),
            ofile2,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
        (true, false) => paired::reader_seq_refine_paired_read(
            MultiGzDecoder::new(reader1),
            ofile1,
            reader2,
            ofile2,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
        (false, true) => paired::reader_seq_refine_paired_read(
            reader1,
            ofile1,
            MultiGzDecoder::new(reader2),
            ofile2,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
        (false, false) => paired::reader_seq_refine_paired_read(
            reader1,
            ofile1,
            reader2,
            ofile2,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
    }
}
