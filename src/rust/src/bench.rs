use std::fs::File;
use std::io::Cursor;
use std::path::Path;

use anyhow::Result;
use extendr_api::prelude::*;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressFinish};
use memchr::memrchr;
use memmap2::Mmap;

use crate::reader::bytes::{BytesReader, ProgressBarReader};

#[extendr]
fn bench_read(file: &str, chunk_size: usize, mmap: bool) -> std::result::Result<(), String> {
    read_chunk(file, chunk_size, mmap).map_err(|e| format!("{}", e))
}

fn read_chunk(file: &str, mut chunk_size: usize, mmap: bool) -> Result<()> {
    let path: &Path = file.as_ref();
    let reader = File::open(file)?;
    let size = reader.metadata()?.len();
    let style = crate::progress_style()?;
    let pb = ProgressBar::new(size).with_finish(ProgressFinish::Abandon);
    pb.set_prefix(format!("Reading {}", file));
    pb.set_style(style);
    if mmap {
        let map = unsafe { Mmap::map(&reader) }?;
        crate::mmap_advice(&map)?;

        if crate::gz_compressed(path) {
            let reader = BytesReader::new(MultiGzDecoder::new(ProgressBarReader::new(
                Cursor::new(map),
                pb,
            )));
            let mut reader = crate::kractor::koutput::io::KoutputBytesChunkReader::new(reader);
            while let Some(_) = reader.chunk_reader()? {}
        } else {
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
        }
    } else {
        let reader = ProgressBarReader::new(reader, pb);
        if crate::gz_compressed(path) {
            let reader = BytesReader::new(MultiGzDecoder::new(reader));
            let mut reader = crate::kractor::koutput::io::KoutputBytesChunkReader::new(reader);
            while let Some(_) = reader.chunk_reader()? {}
        } else {
            let reader = BytesReader::new(reader);
            let mut reader = crate::kractor::koutput::io::KoutputBytesChunkReader::new(reader);
            while let Some(_) = reader.chunk_reader()? {}
        }
    }
    Ok(())
}

extendr_module! {
    mod bench;
    fn bench_read;
}
