use extendr_api::prelude::*;
mod koutput;
mod reads;

use aho_corasick::AhoCorasick;
use koutput::KOutputProcessor;
use reads::write_matching_reads;

use crate::chunk::ChunkProcessor;

#[extendr]
#[allow(clippy::too_many_arguments)]
fn kractor(
    koutput: &str,
    patterns: Robj,
    ofile: &str,
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
    let pattern_vec = patterns
        .as_str_vector()
        .ok_or("`patterns` must be a character vector")?;
    let matcher = AhoCorasick::new(pattern_vec).map_err(|e| {
        format!("Failed to create Aho-Corasick automaton: {}", e)
    })?;
    KOutputProcessor::new(matcher).chunk_io(
        koutput,
        ofile,
        read_buffer,
        write_buffer,
        parse_buffer,
        read_queue,
        write_queue,
        threads,
    )?;
    write_matching_reads(ofile, fq1, ofile1, fq2, ofile2, write_buffer)
}

extendr_module! {
    mod kractor;
    fn kractor;
}
