#' Extracting microbiome kraken2 output and sequence reads
#'
#' @inheritParams blit::kraken2
#' @inheritParams extract_kraken_output
#' @inheritParams extract_kraken_reads
#' @export
kraken2_microbiome <- function(fq1, ..., fq2 = NULL, db = NULL,
                               taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                               kreport = "kraken_report.txt",
                               koutput = "kraken_output.txt",
                               classified_out = "classified.fq",
                               unclassified_out = NULL,
                               microbiome_koutput = "kraken_microbiome_output.txt",
                               microbiome_reads = NULL,
                               threads = NULL,
                               odir = getwd(),
                               kraken2 = NULL,
                               envpath = NULL,
                               overwrite = FALSE) {
    assert_string(fq1, allow_empty = FALSE)
    assert_string(fq2, allow_empty = FALSE, allow_null = TRUE)
    taxon <- as.character(taxon)
    if (length(taxon) == 0L) {
        cli::cli_abort("empty {.arg taxon} provided")
    }
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    assert_string(classified_out, allow_empty = FALSE)
    assert_string(unclassified_out, allow_empty = FALSE, allow_null = TRUE)
    assert_string(odir, allow_empty = TRUE, allow_null = TRUE)
    if (is.null(odir) || !nzchar(odir)) odir <- getwd()
    assert_number_whole(threads,
        min = 1, max = parallel::detectCores(),
        allow_null = TRUE
    )
    assert_string(db, allow_empty = FALSE, allow_null = TRUE)
    reads <- c(fq1, fq2)
    if (is.null(microbiome_reads)) {
        if (is_scalar(reads)) {
            microbiome_reads <- "kraken_microbiome_reads.fa"
        } else {
            microbiome_reads <- sprintf(
                "kraken_microbiome_reads_%d.fa", seq_along(reads)
            )
        }
    } else {
        microbiome_reads <- as.character(microbiome_reads)
        if (length(microbiome_reads) != length(reads)) {
            cli::cli_abort(paste(
                "{.arg microbiome_reads} must have the same length of the",
                "sequence file {.arg fq1} and {.arg fq2}"
            ))
        }
    }

    threads <- threads %||% parallel::detectCores()
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
    microbiome_koutput_file <- file.path(odir, microbiome_koutput)
    if (overwrite || all(file.exists(microbiome_koutput_file))) {
        cli::cli_inform("Extracting microbiome kraken output")
        taxids <- extract_taxids(file.path(odir, kreport), taxon = taxon)
        # then we extract the kraken2 output for these taxon.
        extract_kraken_output(
            file.path(odir, koutput),
            taxids = taxids,
            ofile = microbiome_koutput,
            odir = odir
        )
    } else {
        cli::cli_inform("Using file exist: {.path {microbiome_koutput_file}}")
    }
    microbiome_reads_file <- file.path(odir, microbiome_reads)
    if (overwrite || all(file.exists(microbiome_reads_file))) {
        cli::cli_inform("Extracting microbiome sequence reads")
        extract_kraken_reads(
            kraken_out = file.path(odir, microbiome_koutput),
            reads = reads,
            ofile = microbiome_reads,
            odir = odir,
            threads = 10L,
            # try to change `seqkit` argument into your seqkit path
            # If `NULL`, the internal will detect it in your `PATH` environment variable
            seqkit = NULL
        )
    } else {
        cli::cli_inform("Using file{?s} exist: {.path {microbiome_reads_file}}")
    }

    # styler: off
    if (!overwrite && all(file.exists(c(kraken2_files, microbiome_koutput_file,
                                        microbiome_reads_file)))) {
        cli::cli_warn(paste(
            "Nothing to do, please set {.code overwrite = TRUE}",
            "if you want to overwrite the files exist"
        ))
    }
    # styler: on
    cli::cli_inform("Finished!")
}
