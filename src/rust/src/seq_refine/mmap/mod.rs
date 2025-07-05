use std::fs::File;
use std::io::Cursor;
use std::path::Path;

use anyhow::{anyhow, Result};
use flate2::read::MultiGzDecoder;
use indicatif::{MultiProgress, ProgressBar, ProgressFinish};
use memmap2::Mmap;

mod paired;
mod single;

use crate::reader::bytes::ProgressBarReader;
use crate::reader::slice::SliceProgressBarReader;
use crate::seq_action::*;
use crate::seq_refine::io;

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
    let path: &Path = fq1.as_ref();
    let file = File::open(path)?;
    let map = unsafe { Mmap::map(&file) }?;
    crate::mmap_advice(&map)?;

    let style = crate::progress_style()?;
    let pb = ProgressBar::new(map.len() as u64).with_finish(ProgressFinish::Abandon);
    pb.set_prefix("Parsing fq1");
    pb.set_style(style);

    if crate::gz_compressed(path) {
        io::single::reader_seq_refine_single_read(
            ProgressBarReader::new(MultiGzDecoder::new(Cursor::new(map)), pb),
            ofile1,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        )
    } else {
        let mut reader = SliceProgressBarReader::new(&map);
        reader.set_label("fq1");

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
    let path1: &Path = fq1.as_ref();
    let file1 = File::open(path1)?;
    let map1 = unsafe { Mmap::map(&file1) }?;
    crate::mmap_advice(&map1)?;

    let path2: &Path = fq2.as_ref();
    let file2 = File::open(path2)?;
    let map2 = unsafe { Mmap::map(&file2) }?;
    crate::mmap_advice(&map2)?;

    let style = crate::progress_style()?;
    let progress = MultiProgress::new();
    let pb1 =
        progress.add(ProgressBar::new(map1.len() as u64).with_finish(ProgressFinish::Abandon));
    pb1.set_prefix("Parsing fq1");
    pb1.set_style(style.clone());

    let pb2 =
        progress.add(ProgressBar::new(map2.len() as u64).with_finish(ProgressFinish::Abandon));
    pb2.set_prefix("Parsing fq2");
    pb2.set_style(style);
    let actions = SubseqPairedActions::new(actions1, actions2);
    match (crate::gz_compressed(path1), crate::gz_compressed(path2)) {
        (true, true) => io::paired::reader_seq_refine_paired_read(
            ProgressBarReader::new(MultiGzDecoder::new(Cursor::new(map1)), pb1),
            ofile1,
            ProgressBarReader::new(MultiGzDecoder::new(Cursor::new(map2)), pb2),
            ofile2,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
        (true, false) => io::paired::reader_seq_refine_paired_read(
            ProgressBarReader::new(MultiGzDecoder::new(Cursor::new(map1)), pb1),
            ofile1,
            ProgressBarReader::new(Cursor::new(map2), pb2),
            ofile2,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
        (false, true) => io::paired::reader_seq_refine_paired_read(
            ProgressBarReader::new(Cursor::new(map1), pb1),
            ofile1,
            ProgressBarReader::new(MultiGzDecoder::new(Cursor::new(map2)), pb2),
            ofile2,
            actions,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        ),
        (false, false) => {
            let mut reader1 = SliceProgressBarReader::new(&map1);
            reader1.set_label("fq1");

            let mut reader2 = SliceProgressBarReader::new(&map2);
            reader2.set_label("fq2");

            #[cfg(not(test))]
            reader1.attach_bar(pb1);

            #[cfg(not(test))]
            reader2.attach_bar(pb2);

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
    }
}
