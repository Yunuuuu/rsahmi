run_kraken2 <- function(fq1, ..., fq2 = NULL, sample = NULL, out_dir = getwd(),
                        ncbi_blast_path = NULL, kraken2_db = NULL,
                        mpa_report = TRUE, cores = parallel::detectCores(),
                        cmd = NULL, python_cmd = NULL, sys_args = list()) {
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
        "--classified-out",
        file_path(out_dir,
            if (is.null(fq2)) sample else sprintf("%s#", sample),
            ext = "fq"
        ),
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
    if (mpa_report) {
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
    if (mpa_report) {
        mpa_status <- report_mpa(
            file_path(out_dir, sample, ext = "kraken.report.txt"),
            cmd = python_cmd,
            sample = sample, out_dir = out_dir,
            sys_args = sys_args
        )
        status <- c(status, mpa_status)
    }
    status
}

report_mpa <- function(report, cmd = NULL, sample = NULL, out_dir = getwd(), sys_args = list()) {
    sample <- sample %||% sub(
        "\\.kraken\\.report\\.txt$", "", basename(report),
        perl = TRUE
    )
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    run_command(
        args = c(
            "-f1-3,6-8",
            report, ">",
            file_path(out_dir, sample, ext = "kraken.report.std.txt")
        ),
        cmd = NULL,
        name = "cut",
        sys_args = sys_args
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
        mpa_cmd,
        "-r",
        file_path(out_dir, sample, ext = "kraken.report.std.txt"),
        "-o",
        file_path(out_dir, sample, ext = "kraken.report.mpa.txt"),
        "--intermediate-ranks"
    )
    if (nzchar(Sys.which("python2")) == 1L && nzchar(Sys.which("python3")) == 0L) {
        name <- "python2"
    } else {
        name <- "python3"
    }
    run_command(
        args = args,
        cmd = cmd,
        name = name,
        sys_args = sys_args
    )
}
