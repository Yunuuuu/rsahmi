# we must enable the `bench` features to run following code
Sys.setenv(RSAHMI_FEATURES = "bench")
# pak::pak("Yunuuuu/rsahmi")
pak::pak(pkgbuild::build(binary = TRUE))

reader_bench_out_1m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 1 * 1024 * 1024,
        mmap = FALSE,
        pprof = "bench/pprof_kractor_koutput.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)

reader_bench_out_10m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 10 * 1024 * 1024,
        mmap = FALSE,
        pprof = "bench/pprof_kractor_koutput.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)

mmap_bench_out_1m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 1 * 1024 * 1024,
        pprof = "bench/pprof_mmap_kractor_koutput.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)

mmap_bench_out_10m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 10 * 1024 * 1024,
        pprof = "bench/pprof_mmap_kractor_koutput.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)

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

mmap_time_10m <- bench::mark(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        chunk_size = 10 * 1024 * 1024,
        pprof = "bench/pprof_mmap_kractor_koutput.svg"
    ),
    iterations = 1, memory = FALSE, check = FALSE
)
