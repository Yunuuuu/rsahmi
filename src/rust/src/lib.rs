use extendr_api::prelude::*;

mod batchsender;
mod fastq_reader;
mod koutput_reads;
mod kractor;
mod kreport;
mod parser;
mod reader;
mod reader0;
mod seq_action;
mod seq_range;
mod seq_refine;
mod seq_tag;
pub(crate) mod utils;

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
    use koutput_reads;
    use kractor;
}
