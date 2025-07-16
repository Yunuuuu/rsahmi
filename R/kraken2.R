#' Run Kraken2 for Taxonomic Classification
#'
#' This function runs the Kraken2 classifier on one or two FASTQ files (single-
#' or paired-end). It wraps around the `blit::kraken2()` command interface and
#' adds optional support for environment setup (e.g., Conda or custom `PATH`).
#'
#' @inheritParams blit::kraken2
#' @param db Path to the Kraken2 database. You can download prebuilt databases
#'   from <https://benlangmead.github.io/aws-indexes/k2>, or build your own by
#'   following the instructions at
#'   <https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases>.
#' @param kraken2 Optional. Path to the Kraken2 binary if not in the system
#' `PATH`.
#' @param envpath Optional. Additional directory to prepend to the system `PATH`
#'   environment variable.  Useful when `kraken2` is installed in a non-default
#'   location.
#' @param conda Optional. Name of the Conda environment to activate before
#' running Kraken2.
#' @param condaroot A string specifying the path to the conda root prefix. If
#' not provided, the function searches for the root in the following order:
#'   1. the [option] `blit.conda.root`.
#'   2. the [environment variable][Sys.getenv()] `BLIT_CONDA_ROOT`.
#'   3. the root prefix of [`appmamba()`][blit::appmamba].
#' @return None. This function is called for its side effects.
#'   It produces the following output files in `odir` (or the working directory
#'   if `odir` is `NULL`):
#'   - `kreport`: Kraken2 classification report
#'   - `koutput`: Kraken2 raw classification output
#'   - `classified_out`: FASTQ file of classified reads
#'   - `unclassified_out`: FASTQ file of unclassified reads (if specified)
#' @export
kraken2 <- function(reads, ...,
                    db = NULL,
                    kreport = "kraken_report.txt",
                    koutput = "kraken_output.txt",
                    classified_out = "classified.fq",
                    unclassified_out = NULL,
                    kraken2 = NULL, envpath = NULL,
                    conda = NULL, condaroot = NULL, odir = NULL) {
    reads <- as.character(reads)
    if (length(reads) < 1L || length(reads) > 2L) {
        cli::cli_abort("{.arg reads} must be of length 1 or 2")
    }
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    assert_string(classified_out, allow_empty = FALSE)
    assert_string(unclassified_out, allow_empty = FALSE, allow_null = TRUE)
    assert_string(db, allow_empty = FALSE, allow_null = TRUE)
    assert_number_whole(threads,
        min = 1, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    if (!is.null(db)) db <- sprintf("--db %s", db)
    threads <- threads %||% parallel::detectCores()
    command <- blit::kraken2(
        reads = reads,
        ...,
        ofile = koutput,
        report = kreport,
        classified_out = classified_out,
        unclassified_out = unclassified_out,
        db,
        "--report-minimizer-data",
        sprintf("--threads %d", threads),
        odir = odir,
        kraken2 = kraken2
    )
    command <- blit::cmd_envpath(command, envpath)
    command <- blit::cmd_condaenv(command, conda, root = condaroot)
    blit::cmd_run(command, spinner = TRUE, verbose = TRUE)
}
