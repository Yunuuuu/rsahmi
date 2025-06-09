use extendr_api::prelude::*;
mod koutput;
mod reads;

use ahash::AHashSet as HashSet;
use aho_corasick::AhoCorasick;
use koutput::KOutputProcessor;
use reads::{read_sequence_id_from_koutput, ReadsProcessor};

use crate::chunk::ChunkProcessor;

#[extendr]
#[allow(clippy::too_many_arguments)]
fn kractor_koutput(
    koutput: &str,
    patterns: Robj,
    ofile: &str,
    read_buffer: usize,
    write_buffer: usize,
    parse_buffer: usize,
    read_queue: usize,
    write_queue: usize,
    threads: usize,
) -> std::result::Result<(), String> {
    let pattern_vec = patterns
        .as_str_vector()
        .ok_or("`patterns` must be a character vector")?;
    let matcher = AhoCorasick::new(pattern_vec).map_err(|e| {
        format!("Failed to create Aho-Corasick automaton: {}", e)
    })?;
    rprintln!("Extracting matching kraken2 output from: {}", koutput);
    KOutputProcessor::new(matcher).chunk_io(
        koutput,
        ofile,
        read_buffer,
        write_buffer,
        parse_buffer,
        read_queue,
        write_queue,
        threads,
    )
}

#[extendr]
#[allow(clippy::too_many_arguments)]
fn kractor_reads(
    koutput: &str,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    read_buffer: usize,
    write_buffer: usize,
    parse_buffer: usize,
    read_queue: usize,
    write_queue: usize,
    threads: usize,
) -> std::result::Result<(), String> {
    rprintln!("Extracting sequence IDs");
    let ids = read_sequence_id_from_koutput(koutput, 126 * 1024)
        .map_err(|e| format!("Failed to read sequence IDs: {}", e))?;
    let id_sets = ids
        .iter()
        .map(|id| id.as_slice())
        .collect::<HashSet<&[u8]>>();
    rprintln!("Extracting the matching sequence from: {}", fq1);
    // Don't use multiple threads for the reads processing
    // as it will disturb the order of reads
    ReadsProcessor::new(&id_sets).chunk_io(
        fq1,
        ofile1,
        read_buffer,
        write_buffer,
        parse_buffer,
        read_queue,
        write_queue,
        threads,
    )?;
    if let (Some(in_file), Some(out_file)) = (fq2, ofile2) {
        rprintln!("Extracting the matching sequence from: {}", in_file);
        ReadsProcessor::new(&id_sets).chunk_io(
            in_file,
            out_file,
            read_buffer,
            write_buffer,
            parse_buffer,
            read_queue,
            write_queue,
            threads,
        )?;
    }
    Ok(())
}

#[extendr]
#[allow(clippy::too_many_arguments)]
#[cfg(feature = "bench")]
fn pprof_kractor_koutput(
    koutput: &str,
    patterns: Robj,
    ofile: &str,
    read_buffer: usize,
    write_buffer: usize,
    parse_buffer: usize,
    read_queue: usize,
    write_queue: usize,
    threads: usize,
    pprof_file: &str,
) -> std::result::Result<(), String> {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .build()
        .map_err(|e| format!("cannot create profile guard {:?}", e))?;
    let out = kractor_koutput(
        koutput,
        patterns,
        ofile,
        read_buffer,
        write_buffer,
        parse_buffer,
        read_queue,
        write_queue,
        threads,
    );
    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create(pprof_file).unwrap();
        let mut options = pprof::flamegraph::Options::default();
        options.image_width = Some(2500);
        report.flamegraph_with_options(file, &mut options).unwrap();
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
}
