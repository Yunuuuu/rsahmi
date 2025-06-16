# we must enable the `bench` features to run following code
Sys.setenv(RSAHMI_FEATURES = "bench")
pak::pak(pkgbuild::build(binary = TRUE), ask = FALSE)

system.time(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        # read_buffer = 1024L * 1024L,
        # threads = 5L,
        read_queue = Inf,
        # write_queue = 1000,
        pprof = "bench/pprof_kractor_koutput.svg"
    )
)

rsahmi:::rust_bench_kractor_koutput_reader(
    "bench/data/CNP000460_P01N_output.txt",
    
)

system.time(
    rsahmi:::rust_kractor_reads(
        "bench/data/kraken_microbiome_output.txt",
        "bench/data/CNP000460_P01N_read1.fq",
        extract_reads = "bench/data/kraken_microbiome_reads_1.fa",
        pprof = "bench/pprof_kractor_reads.svg"
    )
)
