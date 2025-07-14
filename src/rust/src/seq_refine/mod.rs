use anyhow::{anyhow, Result};
use extendr_api::prelude::*;
use indicatif::{MultiProgress, ProgressBar, ProgressFinish};

mod paired;
mod single;

use crate::seq_action::*;
use crate::utils::*;

#[extendr]
fn seq_refine(
    fq1: &str,
    ofile1: Option<&str>,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    actions1: Robj,
    actions2: Robj,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: i32,
    nqueue: Option<usize>,
    threads: usize,
) -> std::result::Result<(), String> {
    let actions1 = robj_to_seq_actions(&actions1).map_err(|e| format!("actions1: {}", e))?;
    let actions2 = robj_to_seq_actions(&actions2).map_err(|e| format!("actions2: {}", e))?;
    let threads = threads.max(1); // always use at least one thread
    if let Some(fq2) = fq2 {
        seq_refine_paired_read(
            fq1,
            ofile1,
            fq2,
            ofile2,
            actions1,
            actions2,
            batch_size,
            chunk_bytes,
            compression_level,
            nqueue,
            threads,
        )
        .map_err(|e| format!("{}", e))
    } else {
        seq_refine_single_read(
            fq1,
            ofile1,
            actions1,
            batch_size,
            chunk_bytes,
            compression_level,
            nqueue,
            threads,
        )
        .map_err(|e| format!("{}", e))
    }
}

#[extendr]
#[cfg(feature = "bench")]
fn pprof_seq_refine(
    fq1: &str,
    ofile1: Option<&str>,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    actions1: Robj,
    actions2: Robj,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: i32,
    nqueue: Option<usize>,
    threads: usize,
    pprof_file: &str,
) -> std::result::Result<(), String> {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(2000)
        .build()
        .map_err(|e| format!("cannot create profile guard {:?}", e))?;
    let out = seq_refine(
        fq1,
        ofile1,
        fq2,
        ofile2,
        actions1,
        actions2,
        batch_size,
        chunk_bytes,
        compression_level,
        nqueue,
        threads,
    );
    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create(pprof_file)
            .map_err(|e| format!("Failed to create file {}: {}", pprof_file, e))?;
        let mut options = pprof::flamegraph::Options::default();
        options.image_width = Some(2500);
        report
            .flamegraph_with_options(file, &mut options)
            .map_err(|e| format!("Failed to write flamegraph to {}: {}", pprof_file, e))?;
    };
    out
}

fn seq_refine_single_read(
    fq1: &str,
    ofile1: Option<&str>,
    actions: Option<SubseqActions>,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: i32,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let ofile1 = ofile1.ok_or_else(|| anyhow!("No output file specified."))?;
    let actions = actions.ok_or_else(|| anyhow!("No sequence actions were specified."))?;
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

    single::seq_refine_single_read(
        &fq1,
        Some(pb1),
        &ofile1,
        Some(pb2),
        &actions,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    )
}

fn seq_refine_paired_read(
    fq1: &str,
    ofile1: Option<&str>,
    fq2: &str,
    ofile2: Option<&str>,
    actions1: Option<SubseqActions>,
    actions2: Option<SubseqActions>,
    batch_size: usize,
    chunk_bytes: usize,
    compression_level: i32,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    if ofile1.is_none() && ofile2.is_none() {
        return Err(anyhow!("No output file specified."));
    }
    if actions1.is_none() && actions2.is_none() {
        return Err(anyhow!(
            "No sequence actions were specified. Please provide at least one action to proceed"
        ));
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

    let actions = SubseqPairedActions::new(actions1, actions2);
    paired::seq_refine_paired_read(
        fq1,
        Some(pb1),
        fq2,
        Some(pb3),
        ofile1,
        pb2,
        ofile2,
        pb4,
        &actions,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    )
}

#[cfg(not(feature = "bench"))]
extendr_module! {
    mod seq_refine;
    fn seq_refine;
}

#[cfg(feature = "bench")]
extendr_module! {
    mod seq_refine;
    fn seq_refine;
    fn pprof_seq_refine;
}
