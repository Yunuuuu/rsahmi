use std::fs::File;

use anyhow::Result;
use extendr_api::prelude::*;
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
#[cfg(unix)]
use memmap2::Advice;
use memmap2::Mmap;

mod mmap;

use mmap::mmap_seq_refine_single_read;

// use crate::reader::bytes::BytesProgressBarReader;
use crate::reader::slice::SliceProgressBarReader;
use crate::seq_action::*;

#[extendr]
fn seq_refine(
    fq1: &str,
    ofile1: &str,
    _fq2: Option<&str>,
    _ofile2: Option<&str>,
    actions1: Robj,
    _actions2: Option<Robj>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    mmap: bool,
    threads: usize,
) -> std::result::Result<(), String> {
    let rayon_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| format!("Failed to initialize rayon thread pool: {:?}", e))?;
    let actions1 =
        robj_to_seq_range_actions(&actions1, "actions1").map_err(|e| format!("{}", e))?;
    rayon_pool
        .install(|| {
            seq_refine_single_read(
                fq1,
                ofile1,
                actions1,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
                mmap,
            )
        })
        .map_err(|e| format!("{}", e))
}

fn seq_refine_single_read(
    fq1: &str,
    ofile1: &str,
    actions: SubseqActions,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    _mmap: bool,
) -> Result<()> {
    let file = File::open(fq1)?;
    let progress_style = ProgressStyle::with_template(
        "{prefix:.bold.green} {wide_bar:.cyan/blue} {decimal_bytes}/{decimal_total_bytes} [{elapsed_precise}] {decimal_bytes_per_sec} ({eta})",
    )?;

    let map = unsafe { Mmap::map(&file) }?;
    #[cfg(unix)]
    if Advice::WillNeed.is_supported() {
        map.advise(Advice::WillNeed)?;
    } else if Advice::Sequential.is_supported() {
        map.advise(Advice::Sequential)?;
    }
    let mut reader = SliceProgressBarReader::new(&map);
    reader.set_label("reads");
    let pb = ProgressBar::new(map.len() as u64).with_finish(ProgressFinish::Abandon);
    pb.set_prefix("Parsing reads");
    pb.set_style(progress_style);

    #[cfg(not(test))]
    reader.attach_bar(pb);

    mmap_seq_refine_single_read(
        reader,
        ofile1,
        actions,
        chunk_size,
        buffer_size,
        batch_size,
        nqueue,
    )
}

extendr_module! {
    mod seq_refine;
    fn seq_refine;
}
