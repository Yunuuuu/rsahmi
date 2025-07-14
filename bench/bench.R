# we must enable the `bench` features to run following code
Sys.setenv(RSAHMI_FEATURES = "bench") # TO use pprof
# pak::pak("Yunuuuu/rsahmi")
pak::pak(pkgbuild::build(binary = TRUE))
system.time(rsahmi:::rust_seq_refine(
    c(
        "/home/yun/SINLGE_CELL/Rawdata/FL/FL-2_S1_L005_R1_001.fastq.gz",
        "/home/yun/SINLGE_CELL/Rawdata/FL/FL-2_S1_L005_R2_001.fastq.gz"
    ),
    "bench/data/FL-2_R1_umi.fastq.gz",
    "bench/data/FL-2_R2_umi.fastq.gz",
    umi_action1 = rsahmi::seq_range(end = 12),
    barcode_action1 = rsahmi::seq_range(start = 13, end = 15),
    pprof = "pprof_seq_refine.svg",
    threads = 3L
))
bench::mark(
    kr <- rsahmi:::kraken_report("bench/data/CNP000460_P01N_report.txt")
)
tibble::as_tibble(kr)
rsahmi:::bench_read("bench/data/CNP000460_P01N_output.txt")
rsahmi:::bench_read("bench/data/CNP000460_P01N_output.txt", mmap = FALSE)
rsahmi:::bench_read("bench/data/CNP000460_P01N_read1.fq")
rsahmi:::bench_read("bench/data/CNP000460_P01N_read1.fq", mmap = FALSE)

# - The CNP000460 dataset is used only for benchmarking purposes.
# - Actual analysis should verify whether sequencing is
#   single/paired/ubread(10x).
#
# NOTE:
# ──> Memory-mapped I/O (`mmap = TRUE`) appears faster on repeated runs
#     due to OS-level disk caching (page cache), not actual I/O speed.
# ──> To get unbiased results, always restart R before each benchmark.
# ──────────────────────────────────────────────────────────────────────
# Benchmarking "reader" strategy (no mmap)
# ──────────────────────────────────────────────────────────────────────
# This reads the file via standard buffered I/O using the specified chunk size.
# We use `pprof` profiling and save output for further performance analysis.
reader_bench_out_1m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 1 * 1024 * 1024,
        mmap = FALSE,
        pprof = "bench/pprof_reader_kractor_koutput_1m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(reader_bench_out_1m, "bench/reader_bench_out_1m.rds")

# Benchmark with 10MB chunks (non-mmap reader)
reader_bench_out_10m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 10 * 1024 * 1024,
        mmap = FALSE,
        pprof = "bench/pprof_reader_kractor_koutput_10m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(reader_bench_out_10m, "bench/reader_bench_out_10m.rds")

# ──────────────────────────────────────────────────────────────────────
# Benchmarking "mmap" strategy
# ──────────────────────────────────────────────────────────────────────
mmap_bench_out_1m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 1 * 1024 * 1024,
        mmap = TRUE,
        pprof = "bench/pprof_mmap_kractor_koutput_1m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(mmap_bench_out_1m, "bench/mmap_bench_out_1m.rds")

mmap_bench_out_10m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 10 * 1024 * 1024,
        mmap = TRUE,
        pprof = "bench/pprof_mmap_kractor_koutput_10m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(mmap_bench_out_10m, "bench/mmap_bench_out_10m.rds")

dplyr::bind_rows(
    reader_bench_out_1m,
    reader_bench_out_10m,
    mmap_bench_out_1m,
    mmap_bench_out_10m
)
#  expression             min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc
#   <bch:expr>          <bch:> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>
# 1 "rsahmi:::rust_kra… 10.39m 10.39m   0.00160    22.2MB    0.285     1   178
# 2 "rsahmi:::rust_kra… 16.22m 16.22m   0.00103     3.1MB    0.166     1   162
# 3 "rsahmi:::rust_kra… 12.13m 12.13m   0.00137     3.1MB    0.206     1   150
# 4 "rsahmi:::rust_kra…  2.78m  2.78m   0.00599     3.1MB    0.947     1   158

reader_bench_reads_1m <- bench::mark(
    rsahmi:::rust_kractor_reads(
        "bench/data/kraken_microbiome_output.txt",
        "bench/data/CNP000460_P01N_read1.fq",
        extract_reads = "bench/data/kraken_microbiome_reads.fa",
        chunk_size = 1 * 1024 * 1024,
        mmap = FALSE,
        pprof = "bench/pprof_reader_kractor_reads_1m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(reader_bench_reads_1m, "bench/reader_bench_reads_1m.rds")

reader_bench_reads_10m <- bench::mark(
    rsahmi:::rust_kractor_reads(
        "bench/data/kraken_microbiome_output.txt",
        "bench/data/CNP000460_P01N_read1.fq",
        extract_reads = "bench/data/kraken_microbiome_reads.fa",
        chunk_size = 10 * 1024 * 1024,
        mmap = FALSE,
        pprof = "bench/pprof_reader_kractor_reads_10m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(reader_bench_reads_10m, "bench/reader_bench_reads_10m.rds")

mmap_bench_reads_1m <- bench::mark(
    rsahmi:::rust_kractor_reads(
        "bench/data/kraken_microbiome_output.txt",
        "bench/data/CNP000460_P01N_read1.fq",
        extract_reads = "bench/data/kraken_microbiome_reads.fa",
        chunk_size = 1 * 1024 * 1024,
        mmap = TRUE,
        pprof = "bench/pprof_mmap_kractor_reads_1m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(mmap_bench_reads_1m, "bench/mmap_bench_reads_1m.rds")

mmap_bench_reads_10m <- bench::mark(
    rsahmi:::rust_kractor_reads(
        "bench/data/kraken_microbiome_output.txt",
        "bench/data/CNP000460_P01N_read1.fq",
        extract_reads = "bench/data/kraken_microbiome_reads.fa",
        chunk_size = 10 * 1024 * 1024,
        mmap = TRUE,
        pprof = "bench/pprof_mmap_kractor_reads_10m.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
saveRDS(mmap_bench_reads_10m, "bench/mmap_bench_reads_10m.rds")

system.time(rsahmi:::rust_kractor_reads(
    "bench/data/classified_kraken_microbiome_output.txt",
    c(
        "bench/data/classified_1.fq",
        "bench/data/classified_1.fq"
    ),
    extract_reads = c(
        "bench/data/classified_microbiome_reads_1.fa",
        "bench/data/classified_microbiome_reads_2.fa"
    ),
    chunk_size = 10 * 1024 * 1024,
    mmap = FALSE,
    pprof = "bench/pprof_reader_paired_kractor_reads_10m.svg"
))
