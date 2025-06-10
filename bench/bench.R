# we must enable the `bench` features to run following code
Sys.setenv(RSAHMI_FEATURES = "bench")
pkgbuild::build(binary = TRUE)
pak::pak(
    "local::../rsahmi_0.0.2.9000_R_x86_64-pc-linux-gnu.tar.gz",
    ask = FALSE
)

system.time(
    rsahmi:::rust_kractor_koutput(
        "bench/data/CNP000460_P01N_report.txt",
        "bench/data/CNP000460_P01N_output.txt",
        extract_koutput = "bench/data/kraken_microbiome_output.txt",
        pprof = "bench/pprof_kractor_koutput.svg"
    )
)

system.time(
    rsahmi:::rust_kractor_reads(
        "bench/data/kraken_microbiome_output.txt",
        "bench/data/CNP000460_P01N_read1.fq",
        extract_reads = "bench/data/kraken_microbiome_reads_1.fa",
        pprof = "bench/pprof_kractor_reads.svg"
    )
)
