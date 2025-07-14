use extendr_api::prelude::*;

use crate::utils::*;

pub(crate) mod koutput;
pub(crate) mod reads;

#[extendr]
#[allow(clippy::too_many_arguments)]
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
    let taxonomy = robj_to_option_str(&taxonomy).map_err(|e| format!("'taxonomy' {}", e))?;
    let ranks = robj_to_option_str(&ranks).map_err(|e| format!("'ranks' {}", e))?;
    let taxa = robj_to_option_str(&taxa).map_err(|e| format!("'taxa' {}", e))?;
    let taxids = robj_to_option_str(&taxids).map_err(|e| format!("'taxids' {}", e))?;
    let exclude = robj_to_option_str(&exclude).map_err(|e| format!("'exclude' {}", e))?;
    rayon_pool
        .install(|| {
            koutput::kractor_koutput(
                kreport,
                koutput,
                taxonomy,
                ranks,
                taxa,
                taxids,
                exclude,
                descendants,
                ofile,
                chunk_size,
                buffer_size,
                batch_size,
                nqueue,
                mmap,
            )
        })
        .map_err(|e| format!("{}", e))
}

#[extendr]
#[allow(clippy::too_many_arguments)]
fn kractor_reads(
    koutput: &str,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
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
    rayon_pool.install(|| {
        reads::kractor_reads(
            koutput,
            fq1,
            ofile1,
            fq2,
            ofile2,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
            mmap,
        )
        .map_err(|e| format!("{}", e))
    })
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
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    mmap: bool,
    threads: usize,
    pprof_file: &str,
) -> std::result::Result<(), String> {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(2000)
        .build()
        .map_err(|e| format!("cannot create profile guard {:?}", e))?;
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
        chunk_size,
        buffer_size,
        batch_size,
        nqueue,
        mmap,
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

#[extendr]
#[allow(clippy::too_many_arguments)]
#[cfg(feature = "bench")]
fn pprof_kractor_reads(
    koutput: &str,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    mmap: bool,
    threads: usize,
    pprof_file: &str,
) -> std::result::Result<(), String> {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(2000)
        .build()
        .map_err(|e| format!("cannot create profile guard {:?}", e))?;
    let out = kractor_reads(
        koutput,
        fq1,
        ofile1,
        fq2,
        ofile2,
        chunk_size,
        buffer_size,
        batch_size,
        nqueue,
        mmap,
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

#[cfg(all(test, feature = "bench"))]
mod bench {
    use std::fs::{remove_file, File};
    use std::io::{BufWriter, Read, Write};
    use std::time::Instant;

    use rand::rngs::StdRng;
    use rand::{RngCore, SeedableRng};

    const FILE_PATH: &str = "1gb_test_file.dat";
    const TOTAL_SIZE: usize = 1024 * 1024 * 1024; // 1 GB
    const CHUNK_SIZE: usize = 1024 * 1024; // 1 MB

    #[test]
    fn bench_io() -> std::io::Result<()> {
        let mut rng = StdRng::seed_from_u64(42);
        let mut file = BufWriter::with_capacity(CHUNK_SIZE, File::create(FILE_PATH)?);

        // === Pre-generate all chunks in memory ===
        let mut chunks: Vec<Vec<u8>> = Vec::with_capacity(TOTAL_SIZE / CHUNK_SIZE);
        for _ in 0 .. (TOTAL_SIZE / CHUNK_SIZE) {
            let mut buffer = vec![0u8; CHUNK_SIZE];
            rng.fill_bytes(&mut buffer);
            chunks.push(buffer);
        }

        // === Benchmark only write time ===
        let start = Instant::now();
        for chunk in &chunks {
            file.write_all(chunk)?;
        }
        file.flush()?;
        let elapsed = start.elapsed();

        println!(
            "Write completed: {:.2} MB in {:.2?} ({:.2} MB/s | {:.2} GB/s)",
            TOTAL_SIZE as f64 / 1024.0 / 1024.0,
            elapsed,
            TOTAL_SIZE as f64 / 1024.0 / 1024.0 / elapsed.as_secs_f64(),
            TOTAL_SIZE as f64 / 1024.0 / 1024.0 / 1024.0 / elapsed.as_secs_f64()
        );

        // === Benchmark only read time ===
        let mut file = File::open(FILE_PATH)?;
        // let mut reader = BufReader::with_capacity(CHUNK_SIZE, file);
        let mut buffer = vec![0u8; CHUNK_SIZE];

        let mut total_read = 0;
        let start = Instant::now();

        loop {
            let bytes_read = file.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }
            total_read += bytes_read;
        }

        let elapsed = start.elapsed();
        println!(
            "Read completed: {:.2} MB in {:.2?} ({:.2} MB/s | {:.2} GB/s)",
            total_read as f64 / 1024.0 / 1024.0,
            elapsed,
            total_read as f64 / 1024.0 / 1024.0 / elapsed.as_secs_f64(),
            total_read as f64 / 1024.0 / 1024.0 / 1024.0 / elapsed.as_secs_f64()
        );

        assert_eq!(total_read, TOTAL_SIZE);
        remove_file(FILE_PATH)?;

        Ok(())
    }
}
