Sys.setenv(RUST_FEATURES = "bench")
pkgbuild::build(binary = TRUE)
pak::pak(
    "local::../rsahmi_0.0.2.9000_R_x86_64-pc-linux-gnu.tar.gz",
    ask = FALSE
)
rsahmi:::rust_kractor_koutput(
    "bench/data/CNP000460_P01N_report.txt",
    "bench/data/CNP000460_P01N_output.txt"
)
