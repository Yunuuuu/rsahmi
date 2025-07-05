use extendr_api::prelude::*;

mod io;
mod mmap;

// use crate::reader::bytes::BytesProgressBarReader;
use crate::seq_action::*;

#[extendr]
fn seq_refine(
    fq1: &str,
    ofile1: Option<&str>,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    actions1: Robj,
    actions2: Robj,
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
    let actions1 = robj_to_seq_actions(&actions1, "actions1").map_err(|e| format!("{}", e))?;
    let actions2 = robj_to_seq_actions(&actions2, "actions2").map_err(|e| format!("{}", e))?;

    rayon_pool
        .install(|| {
            if mmap {
                mmap::mmap_seq_refine(
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
                io::reader_seq_refine(
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
            }
        })
        .map_err(|e| format!("{}", e))
}

extendr_module! {
    mod seq_refine;
    fn seq_refine;
}
