use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use extendr_api::prelude::*;
use indicatif::style::TemplateError;
use indicatif::ProgressStyle;
#[cfg(unix)]
use memmap2::Advice;
use memmap2::Mmap;

mod batchsender;
mod bench;
// mod koutput_reads;
mod fastq_reader;
mod kractor;
mod kreport;
mod parser;
mod reader;
mod reader0;
mod seq_action;
mod seq_refine;

pub(crate) const BLOCK_SIZE: usize = 4 * 1024 * 1024;

pub(crate) fn new_channel<T>(nqueue: Option<usize>) -> (Sender<T>, Receiver<T>) {
    if let Some(queue) = nqueue {
        bounded(queue)
    } else {
        unbounded()
    }
}

pub(crate) fn mmap_advice(map: &Mmap) -> anyhow::Result<()> {
    #[cfg(unix)]
    if Advice::WillNeed.is_supported() {
        map.advise(Advice::WillNeed)?;
    } else if Advice::Sequential.is_supported() {
        map.advise(Advice::Sequential)?;
    }
    Ok(())
}

pub(crate) fn progress_reader_style() -> std::result::Result<ProgressStyle, TemplateError> {
    ProgressStyle::with_template(
        "{prefix:.bold.cyan/blue} {decimal_bytes}/{decimal_total_bytes} {spinner:.green} [{elapsed_precise}] {decimal_bytes_per_sec} (ETA {eta})",
    )
}

pub(crate) fn progress_writer_style() -> std::result::Result<ProgressStyle, TemplateError> {
    ProgressStyle::with_template(
        "{prefix:.bold.cyan/blue} {decimal_bytes} {spinner:.green} {decimal_bytes_per_sec}",
    )
}

pub(crate) fn gz_compressed(path: &std::path::Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .map_or(false, |s| s.eq_ignore_ascii_case("gz"))
}

// https://extendr.github.io/extendr/extendr_api/#returning-resultt-e-to-r
// https://github.com/extendr/extendr/blob/master/extendr-api/src/robj/into_robj.rs#L100
// The memory-safe way to do error handling with extendr is to return a Result<T, E> to R.
// By default, any Err will trigger a panic! on the rust side which unwinds the stack.
// The rust error trace will be printed to stderr, not R terminal. Any Ok value is returned as is.
//
// 1. we must wrap (Wrap) all third-part object, only local struc or enum will be transformed into R object
// 2. Define a error object, and implement a from method for `Robj`.
// 3. extendr don't automatically transform Struct (Or enum) object
//    in `function` into Robj but only transform in `impl`.

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
// For methods, we'll call it directly with R function `call_rust_method`
extendr_module! {
    mod rsahmi;
    use kreport;
    use seq_refine;
    use kractor;
    use bench;
}
