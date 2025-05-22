#' Extracting microbiome kraken2 output and sequence reads
#'
#' @inheritParams blit::kraken2
#' @inheritParams kractor
#' @export
kraken_taxon <- function(fq1, ..., fq2 = NULL, db = NULL,
                         kreport = "kraken_report.txt",
                         koutput = "kraken_output.txt",
                         classified_out = "classified.fq",
                         unclassified_out = NULL,
                         extract_koutput = NULL,
                         extract_reads = NULL,
                         taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                         buffer_size = NULL, threads = NULL,
                         odir = getwd(),
                         kraken2 = NULL, envpath = NULL,
                         overwrite = FALSE) {
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
    assert_number_whole(buffer_size, min = 1, allow_null = TRUE)
    assert_number_whole(threads,
        min = 1, max = parallel::detectCores(),
        allow_null = TRUE
    )
    threads <- threads %||% parallel::detectCores()
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
            buffer_size = buffer_size,
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
