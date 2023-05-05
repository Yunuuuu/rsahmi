#' Running Kraken2
#' @description 
#' Metagenomic classification of paired-end reads from single-cell RNA
#' sequencing fastq files can be performed using any k-mer based mapper that
#' identifies a taxonomic ID for each k-mer and read. However, `SAHMI` is
#' optimized to run with `Kraken2Uniq`, which finds exact matches of candidate
#' 35-mer genomic substrings to the lowest common ancestor of genomes in a
#' reference metagenomic database. It is essential that all realistically
#' possible genomes are included as mapping references at this stage (e.g. host,
#' known vectors, etc.), or that host mappable reads are excluded. The required
#' outputs from this step are: a Kraken summary report with sample level
#' metagenomic counts, a Kraken output file with read and k-mer level taxonomic
#' classifications, an MPA-style report, and raw sequencing fastq files with
#' taxonomic classification for each read. Please see
#' <https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown> for
#' more details on installation and usage of Kraken2/KrakenUniq. 
#' @param fq1,fq2 Path to fastq 1 file.
#' @param sample A string, sample name, will be used to create output file name.
#' @param out_dir The path of output directory.
#' @param ncbi_blast_path The path of ncbi-blastx.
#' @param kraken2_db The path of kraken database.
#' @param mpa_report A scalar logical indicates whethter reporting mpa-style
#'   results. 
#' @param cores The number of cores to use.
#' @param cmd The path of `kraken2` command.
#' @param kraken2_args Other arguments passed to `kraken2`.
#' @param python_cmd The path of `python2` or `python3`.
#' @param sys_args Other arguments passed to [system2].
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @export 
run_kraken2 <- function(fq1, fq2 = NULL, sample = NULL, out_dir = getwd(),
                        ncbi_blast_path = NULL, kraken2_db = NULL,
                        mpa_report = TRUE, cores = availableCores(),
                        cmd = NULL, kraken2_args = character(),
                        python_cmd = NULL, sys_args = list()) {
    sample <- sample %||% sub("_0*[12]?\\.(fastq|fq)(\\.gz)?$", "",
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
    # https://github.com/DerrickWood/kraken2/wiki/Manual
    # Usage of --paired also affects the --classified-out and
    # --unclassified-out options; users should provide a # character in the
    # filenames provided to those options, which will be replaced by kraken2
    # with "_1" and "_2" with mates spread across the two files
    # appropriately. For example:
    if (is.null(fq2)) {
        classified_out_file <- sample
    } else {
        classified_out_file <- sprintf("%s#", sample)
    }
    args <- c(
        args,
        "--classified-out",
        file_path(out_dir, classified_out_file, ext = "fq"),
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
    args <- c(args, kraken2_args, fq1, fq2)
    if (mpa_report) {
        kraken2_sys_args <- sys_args
        kraken2_sys_args$wait <- TRUE
    }
    if (!is.null(ncbi_blast_path)) {
        cli::cli_alert("Adding {.field blast} into {.var PATH}")
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
    if (nzchar(Sys.which("python2")) && !nzchar(Sys.which("python3"))) {
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
