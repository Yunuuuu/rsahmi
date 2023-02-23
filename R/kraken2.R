run_kraken2 <- function(
    fq1, fq2 = NULL, sample = NULL, out_dir = getwd(),
    ncbi_blast_path = NULL, kraken2_db = NULL,
    report_mpa = TRUE, cores = parallel::detectCores(),
    ..., cmd = NULL, sys_args = list()) {
    sample <- sample %||% sub("_1[.fastq]*$", "", basename(fq1), perl = TRUE)
    args <- handle_arg(cores, "--threads", "%d")
    args <- c(args, handle_arg(kraken2_db, "--db", "%s"))
    args <- c(args, handle_arg(!is.null(fq2), "--paired"))
    args <- c(args, "--use-names", "--report-minimizer-data")
    args <- c(
        args,
        "--classified-out", normalizePath(file_path(out_dir, sample)), "#.fq",
        "--output",
        normalizePath(file_path(out_dir, sample, ext = "kraken.output.txt")),
        "--report",
        normalizePath(file_path(out_dir, sample, ext = "kraken.report.txt"))
    )
    args <- c(args, ..., fq1, fq2)
    if (report_mpa) {
        kraken2_sys_args <- sys_args
        kraken2_sys_args$wait <- TRUE
    }
    if (!is.null(ncbi_blast_path)) {
        old_path <- Sys.getenv("PATH")
        on.exit(Sys.setenv(PATH = old_path))
        Sys.setenv(PATH = paste(old_path,
            normalizePath(ncbi_blast_path),
            sep = ":"
        ))
    }
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }
    status <- run_command(
        args,
        cmd = cmd, name = "kraken2",
        sys_args = kraken2_sys_args
    )
    if (report_mpa) {
        mpa_status <- report_mpa(
            file_path(out_dir, sample, ext = "kraken.report.txt"),
            sample = sample, out_dir = out_dir,
            sys_args = sys_args
        )
        status <- c(status, mpa_status)
    }
    status
}

report_mpa <- function(report, sample = NULL, out_dir = getwd(), sys_args = list()) {
    sample <- sample %||% sub(
        "\\.kraken\\.report\\.txt$", "", basename(report),
        perl = TRUE
    )
    data <- data.table::fread(report, sep = "\t", header = FALSE)
    data <- data[, c(1:3, 6:8)]
    data.table::fwrite(
        data,
        file = file_path(out_dir, sample, ext = "kraken.report.std.txt"),
        sep = "\t", col.names = FALSE
    )
    mpa_cmd <- system.file(
        "extdata", "kreport2mpa.py",
        package = "rsahmi", mustWork = TRUE
    )
    args <- c(
        "-r",
        normalizePath(
            file_path(out_dir, sample, ext = "kraken.report.std.txt")
        ),
        "-o",
        normalizePath(
            file_path(out_dir, sample, ext = "kraken.report.mpa.txt")
        ),
        "--intermediate-ranks"
    )
    run_command(
        args = args,
        cmd = mpa_cmd,
        sys_args = sys_args
    )
}
