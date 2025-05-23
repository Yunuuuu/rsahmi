use extendr_api::prelude::*;
mod output;
mod reads;

use aho_corasick::AhoCorasick;
use output::write_matching_output;
use reads::write_matching_reads;

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
    io_buffer: usize,
    buffersize: usize,
    batchsize: usize,
    queue_capacity: usize,
    threads: usize,
) -> std::result::Result<(), String> {
    let pattern_vec = patterns
        .as_str_vector()
        .ok_or("`patterns` must be a character vector")?;
    let matcher = AhoCorasick::new(pattern_vec).map_err(|e| {
        format!("Failed to create Aho-Corasick automaton: {}", e)
    })?;
    write_matching_output(
        koutput,
        &matcher,
        ofile,
        io_buffer,
        buffersize,
        batchsize,
        queue_capacity,
        threads,
    )?;
    write_matching_reads(ofile, fq1, ofile1, fq2, ofile2, io_buffer)
}

extendr_module! {
    mod kractor;
    fn kractor;
}
