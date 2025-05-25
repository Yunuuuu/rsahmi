#' Extract microbiome Kraken2 output and sequence reads
#'
#' @inheritParams blit::kraken2
#' @param db Path to the Kraken2 database. You can download prebuilt databases
#'   from <https://benlangmead.github.io/aws-indexes/k2>, or build your own by
#'   following the instructions at
#'   <https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases>.
#' @inheritParams kractor
#' @param envpath String. Additional path to include in the `PATH` environment
#'   variable. Used to locate the `kraken2` executable.
#' @param overwrite Logical. Whether to overwrite existing files in `odir`.
#' @seealso [`kractor()`]
#' @return None. This function generates the following files:
#' - `kreport`: the `kraken2` report file.
#' - `koutput`: the `kraken2` output file.
#' - `classified_out`: the classified sequence file(s) from `kraken2`.
#' - `unclassified_out`: the unclassified sequence file(s) from `kraken2`.
#' - `extract_koutput`: Kraken2 output entries corresponding to the specified
#'   `taxon`, extracted from koutput.
#' - `extract_reads`: Sequence file(s) containing reads assigned to the
#'   specified `taxon`.
#' @export
kraken_taxon <- function(fq1, ..., fq2 = NULL, db = NULL,
                         kreport = "kraken_report.txt",
                         koutput = "kraken_output.txt",
                         classified_out = "classified.fq",
                         unclassified_out = NULL,
                         extract_koutput = NULL,
                         extract_reads = NULL,
                         taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                         read_buffer = NULL, write_buffer = NULL,
                         parse_buffer = NULL,
                         read_queue = NULL, write_queue = NULL,
                         threads = NULL, kraken2 = NULL, envpath = NULL,
                         overwrite = FALSE, odir = getwd()) {
    assert_string(fq1, allow_empty = FALSE)
    assert_string(fq2, allow_empty = FALSE, allow_null = TRUE)
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    assert_string(classified_out, allow_empty = FALSE)
    assert_string(unclassified_out, allow_empty = FALSE, allow_null = TRUE)
    assert_string(db, allow_empty = FALSE, allow_null = TRUE)

    reads <- c(fq1, fq2)
    assert_string(extract_koutput, allow_empty = FALSE, allow_null = TRUE)
    if (is.null(extract_reads)) {
        if (is_scalar(reads)) {
            extract_reads <- "kraken_microbiome_reads.fa"
        } else {
            extract_reads <- sprintf(
                "kraken_microbiome_reads_%d.fa", seq_along(reads)
            )
        }
    } else {
        extract_reads <- as.character(extract_reads)
        if (length(extract_reads) != length(reads)) {
            cli::cli_abort(paste(
                "{.arg extract_reads} must have the same length of the",
                "sequence file {.arg fq1} and {.arg fq2}"
            ))
        }
    }
    taxon <- as.character(taxon)
    if (length(taxon) == 0L) {
        cli::cli_abort("empty {.arg taxon} provided")
    }
    assert_number_whole(read_buffer, min = 1, allow_null = TRUE)
    assert_number_whole(write_buffer, min = 1, allow_null = TRUE)
    assert_number_whole(parse_buffer, min = 1, allow_null = TRUE)
    assert_number_whole(read_queue, min = 1, allow_null = TRUE)
    assert_number_whole(write_queue, min = 1, allow_null = TRUE)
    assert_number_whole(threads,
        min = 1, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    threads <- threads %||% parallel::detectCores()
    assert_bool(overwrite)
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    if (is.null(odir)) odir <- getwd()
    if (!is.null(db)) db <- sprintf("--db %s", db)
    kraken2_files <- file.path(
        odir,
        c(koutput, kreport, classified_out, unclassified_out)
    )
    if (overwrite || !all(file.exists(kraken2_files))) {
        command <- blit::kraken2(
            fq1 = fq1,
            fq2 = fq2,
            ...,
            ofile = koutput,
            report = kreport,
            classified_out = classified_out,
            unclassified_out = unclassified_out,
            "--use-names",
            "--report-minimizer-data",
            sprintf("--threads %d", threads),
            db,
            odir = odir,
            kraken2 = kraken2
        )
        command <- blit::cmd_envpath(command, envpath)
        blit::cmd_run(command, spinner = TRUE, verbose = TRUE)
    } else {
        cli::cli_inform("Using files exist: {.path {kraken2_files}}")
    }
    kractor_files <- file.path(odir, c(extract_koutput, extract_reads))
    if (overwrite || !all(file.exists(kractor_files))) {
        kractor(
            kreport = file.path(odir, kreport),
            koutput = file.path(odir, koutput),
            reads = reads,
            extract_koutput = extract_koutput,
            extract_reads = extract_reads,
            taxon = taxon,
            read_buffer = read_buffer,
            write_buffer = write_buffer,
            parse_buffer = parse_buffer,
            read_queue = read_queue,
            write_queue = write_queue,
            threads = threads,
            odir = odir
        )
    } else {
        cli::cli_inform("Using file exist: {.path {extract_kractor_files}}")
    }

    if (!overwrite && all(file.exists(c(kraken2_files, kractor_files)))) {
        cli::cli_warn(paste(
            "Nothing to do, please set {.code overwrite = TRUE}",
            "if you want to overwrite the files exist"
        ))
    }
    invisible(NULL)
}
