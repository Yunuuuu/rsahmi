run_kraken2 <- function(fq1, ..., fq2 = NULL, sample = NULL, out_dir = getwd(),
                        ncbi_blast_path = NULL, kraken2_db = NULL,
                        report_mpa = TRUE, cores = parallel::detectCores(),
                        cmd = NULL, sys_args = list()) {
    sample <- sample %||% sub("_0*1\\.(fastq|fq)(\\.gz)?$", "",
        basename(fq1),
        perl = TRUE
    )
    args <- handle_arg(cores, "--threads", "%d")
    args <- c(args, handle_arg(kraken2_db, "--db", "%s"))
    args <- c(args, handle_arg(!is.null(fq2), "--paired"))
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    fq1 <- normalizePath(fq1, mustWork = TRUE)
    if (!is.null(fq2)) fq2 <- normalizePath(fq2, mustWork = TRUE)
    out_dir <- normalizePath(out_dir, mustWork = TRUE)
    args <- c(
        args,
        "--classified-out", file_path(out_dir, sample),
        "--output",
        file_path(out_dir, sample, ext = "kraken.output.txt"),
        "--report",
        file_path(out_dir, sample, ext = "kraken.report.txt"),
        "--use-names", "--report-minimizer-data",
        handle_arg(
            grepl("\\.gz$", fq1, perl = TRUE),
            name = "--gzip-compressed"
        ),
        handle_arg(
            grepl("\\.bz2$", fq1, perl = TRUE),
            name = "--bzip2-compressed"
        )
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
            normalizePath(ncbi_blast_path, mustWork = TRUE),
            sep = ":"
        ))
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
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    data.table::fwrite(
        data,
        file = file_path(out_dir, sample, ext = "kraken.report.std.txt"),
        sep = "\t", col.names = FALSE
    )
    mpa_cmd <- system.file(
        "extdata", "kreport2mpa.py",
        package = "rsahmi", mustWork = TRUE
    )
    if (file.access(mpa_cmd, mode = 1L) != 0L) {
        Sys.chmod(mpa_cmd, "555")
    }
    out_dir <- normalizePath(out_dir, mustWork = TRUE)
    args <- c(
        "-r",
        file_path(out_dir, sample, ext = "kraken.report.std.txt"),
        "-o",
        file_path(out_dir, sample, ext = "kraken.report.mpa.txt"),
        "--intermediate-ranks"
    )
    run_command(
        args = args,
        cmd = mpa_cmd,
        sys_args = sys_args
    )
}
