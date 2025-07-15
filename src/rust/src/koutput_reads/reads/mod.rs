use anyhow::Result;
use bytes::Bytes;
use indicatif::{MultiProgress, ProgressBar, ProgressFinish, ProgressStyle};
use libdeflater::CompressionLvl;
use rustc_hash::FxHashMap as HashMap;

use crate::seq_tag::*;
use crate::utils::*;

mod paired;
mod single;
mod stream;

pub(super) fn parse_reads(
    koutmap: &HashMap<Bytes, (Bytes, Bytes, Bytes)>,
    fq1: &str,
    fq2: Option<&str>,
    ofile: &str,
    tag_ranges1: Option<TagRanges>,
    tag_ranges2: Option<TagRanges>,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: CompressionLvl,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let reader_style = progress_reader_style()?;
    let progress = MultiProgress::new();
    let reader_pb1 = progress.add(
        ProgressBar::new(std::fs::metadata(fq1)?.len() as u64).with_finish(ProgressFinish::Abandon),
    );
    reader_pb1.set_prefix("Reading fq1");
    reader_pb1.set_style(reader_style.clone());

    let matching_pb = ProgressBar::no_length().with_finish(ProgressFinish::Abandon);
    matching_pb.set_style(ProgressStyle::with_template(
        "Matching {human_len.bold.cyan/blue} reads {spinner:.green}",
    )?);

    let threads = threads.max(1); // always use at least one thread
    if let Some(fq2) = fq2 {
        let reader_pb2 = progress.add(
            ProgressBar::new(std::fs::metadata(fq2)?.len() as u64)
                .with_finish(ProgressFinish::Abandon),
        );
        reader_pb2.set_prefix("Reading fq2");
        reader_pb2.set_style(reader_style);
        let matching_pb = progress.add(matching_pb);
        paired::parse_paired_read(
            koutmap,
            fq1,
            Some(reader_pb1),
            fq2,
            Some(reader_pb2),
            ofile,
            Some(matching_pb),
            &tag_ranges1,
            &tag_ranges2,
            batch_size,
            chunk_bytes,
            compression_level,
            nqueue,
            threads,
        )
    } else {
        let matching_pb = progress.add(matching_pb);
        single::parse_single_read(
            koutmap,
            fq1,
            Some(reader_pb1),
            ofile,
            Some(matching_pb),
            &tag_ranges1,
            batch_size,
            chunk_bytes,
            compression_level,
            nqueue,
            threads,
        )
    }
}
