use std::fs::File;
use std::io::Read;

use anyhow::Result;
use bytes::BytesMut;
use extendr_api::prelude::*;
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use memchr::memrchr;
use memmap2::Mmap;

#[extendr]
fn bench_read(file: &str, chunk_size: usize, mmap: bool) -> std::result::Result<(), String> {
    read_chunk(file, chunk_size, mmap).map_err(|e| format!("{}", e))
}

fn read_chunk(file: &str, mut chunk_size: usize, mmap: bool) -> Result<()> {
    let mut reader = File::open(file)?;
    let size = reader.metadata()?.len();
    let progress_style = ProgressStyle::with_template(
        "{prefix:.bold.green} {wide_bar:.cyan/blue} {decimal_bytes}/{decimal_total_bytes} [{elapsed_precise}] {decimal_bytes_per_sec} ({eta})",
    )?;
    let pb = ProgressBar::new(size).with_finish(ProgressFinish::Abandon);
    pb.set_prefix(format!("Reading {}", file));
    pb.set_style(progress_style);
    if mmap {
        let map = unsafe { Mmap::map(&reader) }?;
        crate::mmap_advice(&map)?;

        let mut pos = 0;
        while pos < map.len() {
            if pos + chunk_size < map.len() {
                let chunk = &map[pos .. pos + chunk_size];
                if let Some(offset) = memrchr(b'\n', chunk) {
                    let _ = &chunk[..= offset];
                    pb.inc(offset as u64);
                    pos += offset + 1;
                } else {
                    chunk_size *= 2;
                    continue;
                }
            } else {
                let chunk = &map[pos ..];
                pb.inc(chunk.len() as u64);
                pos = map.len();
            }
        }
    } else {
        let mut leftover = BytesMut::new();
        loop {
            let leftover_len = leftover.len();
            let mut buf = BytesMut::with_capacity(chunk_size + leftover.len());
            buf.extend_from_slice(&leftover);
            unsafe { buf.set_len(leftover_len + chunk_size) };
            let nbytes = reader.read(&mut buf[leftover_len ..])?;
            unsafe { buf.set_len(leftover_len + nbytes) };
            if nbytes == 0 {
                if buf.is_empty() {
                    return Ok(());
                } else {
                    pb.inc(buf.len() as u64);
                    leftover = BytesMut::new();
                    continue;
                }
            }
            if let Some(pos) = memrchr(b'\n', &buf) {
                // split at newline
                let _ = buf.split_to(pos + 1);
                pb.inc(pos as u64);
                leftover = buf; // remainder for next round
                continue;
            } else {
                // no newline found — either final chunk or partial
                // leave all in leftover and try again with larger chunk
                chunk_size *= 2;
                leftover = buf;
                continue;
            }
        }
    }
    Ok(())
}

extendr_module! {
    mod bench;
    fn bench_read;
}
