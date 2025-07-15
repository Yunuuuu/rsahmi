#[cfg(feature = "bench")]
use anyhow::Context;
use extendr_api::prelude::*;

mod koutput;
pub(crate) mod reads;

#[extendr]
fn kractor_koutput(
    kreport: &str,
    koutput: &str,
    taxonomy: Robj,
    ranks: Robj,
    taxa: Robj,
    taxids: Robj,
    exclude: Robj,
    descendants: bool,
    ofile: &str,
    compression_level: i32,
    batch_size: usize,
    chunk_bytes: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> std::result::Result<(), String> {
    koutput::kractor_koutput(
        kreport,
        koutput,
        ofile,
        taxonomy,
        ranks,
        taxa,
        taxids,
        exclude,
        descendants,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    )
    .map_err(|e| format!("{:?}", e))
}

#[extendr]
fn kractor_reads(
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
) -> std::result::Result<(), String> {
    reads::kractor_reads(
        koutput,
        fq1,
        ofile1,
        fq2,
        ofile2,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    )
    .map_err(|e| format!("{}", e))
}

#[extendr]
#[cfg(feature = "bench")]
fn pprof_kractor_koutput(
    kreport: &str,
    koutput: &str,
    taxonomy: Robj,
    ranks: Robj,
    taxa: Robj,
    taxids: Robj,
    exclude: Robj,
    descendants: bool,
    ofile: &str,
    compression_level: i32,
    batch_size: usize,
    chunk_bytes: usize,
    nqueue: Option<usize>,
    threads: usize,
    pprof_file: &str,
) -> std::result::Result<(), String> {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(2000)
        .build()
        .with_context(|| format!("cannot create profile guard"))
        .map_err(|e| format!("{:?}", e))?;
    let out = kractor_koutput(
        kreport,
        koutput,
        taxonomy,
        ranks,
        taxa,
        taxids,
        exclude,
        descendants,
        ofile,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    );
    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create(pprof_file)
            .with_context(|| format!("Failed to create file {}", pprof_file))
            .map_err(|e| format!("{:?}", e))?;
        let mut options = pprof::flamegraph::Options::default();
        options.image_width = Some(2500);
        report
            .flamegraph_with_options(file, &mut options)
            .with_context(|| format!("Failed to write flamegraph to {}", pprof_file))
            .map_err(|e| format!("{:?}", e))?;
    };
    out
}

#[extendr]
#[allow(clippy::too_many_arguments)]
#[cfg(feature = "bench")]
fn pprof_kractor_reads(
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
    pprof_file: &str,
) -> std::result::Result<(), String> {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(2000)
        .build()
        .with_context(|| format!("cannot create profile guard"))
        .map_err(|e| format!("{:?}", e))?;
    let out = kractor_reads(
        koutput,
        fq1,
        ofile1,
        fq2,
        ofile2,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    );
    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create(pprof_file)
            .with_context(|| format!("Failed to create file {}", pprof_file))
            .map_err(|e| format!("{:?}", e))?;
        let mut options = pprof::flamegraph::Options::default();
        options.image_width = Some(2500);
        report
            .flamegraph_with_options(file, &mut options)
            .with_context(|| format!("Failed to write flamegraph to {}", pprof_file))
            .map_err(|e| format!("{:?}", e))?;
    };
    out
}

#[cfg(not(feature = "bench"))]
extendr_module! {
    mod kractor;
    fn kractor_koutput;
    fn kractor_reads;
}

#[cfg(feature = "bench")]
extendr_module! {
    mod kractor;
    fn kractor_koutput;
    fn kractor_reads;
    fn pprof_kractor_koutput;
    fn pprof_kractor_reads;
}
