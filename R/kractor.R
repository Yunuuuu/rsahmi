#' Filter Kraken2 Output by Taxon
#'
#' This function filters Kraken2 classification output (`koutput`) by taxonomic
#' identifiers derived from a Kraken2 report (`kreport`). It extracts lines
#' matching the desired `taxonomy`, `ranks`, `taxa`, `taxids`, and `descendants`
#' and writes the filtered results to an output file.
#'
#' @param ofile Optional path to the output file storing the filtered Kraken2
#'   output lines that pass taxonomic and exclusion filters. If `NULL`, a
#'   default filename `"kraken_microbiome_output.txt"` will be used. If the
#'   filename ends with `.gz`, output will be automatically compressed using
#'   gzip.
#' @param taxonomy Character vector. The set of taxonomic groups to include
#' (default: `c("D__Bacteria", "D__Fungi", "D__Viruses")`). This defines the
#' global taxa to consider. If `NULL`, all taxa will be used. If `descendants =
#' TRUE`, only the descendants within these groups will be considered. The
#' selection of taxa can be further refined using the `ranks`, `taxa`, and
#' `taxids` parameters. One of `taxonomy`, `ranks`, `taxa`, or `taxids` must be
#' provided.
#' @param ranks Character vector. The taxonomic ranks to filter by (optional).
#' @param taxa Character vector. Specific taxa to include (optional).
#' @param taxids Character vector. A list of taxid values to filter by
#' (optional).
#' @param exclude A character vector of taxids to exclude sequences from usage.
#' @param descendants Logical. Whether to include descendants of the selected
#' taxa (default: `TRUE`).
#' @inheritParams koutreads
#' @return None. The function generates a filtered Kraken2 output file
#'   containing entries corresponding to the specified `taxonomy`, `ranks`,
#'   `taxa`, `taxids`, and `descendants` extracted from the input `koutput`.
#' @export
kractor_koutput <- function(kreport, koutput, ofile = NULL,
                            taxonomy = c(
                                "D__Bacteria", "D__Fungi", "D__Viruses"
                            ),
                            ranks = NULL,
                            taxa = NULL,
                            taxids = NULL,
                            exclude = NULL,
                            descendants = TRUE,
                            batch_size = NULL, chunk_bytes = NULL,
                            compression_level = 4L,
                            nqueue = NULL, threads = NULL, odir = NULL) {
    rust_kractor_koutput(
        kreport = kreport,
        koutput = koutput,
        ofile = ofile,
        taxonomy = taxonomy,
        ranks = ranks,
        taxa = taxa,
        taxids = taxids,
        exclude = exclude,
        descendants = descendants,
        batch_size = batch_size,
        chunk_bytes = chunk_bytes,
        compression_level = compression_level,
        nqueue = nqueue,
        threads = threads,
        odir = odir
    )
}

#' Extract Reads from Kraken2 Output Based on Classification
#'
#' This function extracts reads corresponding to selected classifications from a
#' Kraken2 output file (`koutput`). Only reads classified to selected taxa will
#' be extracted from the provided sequence file (`reads`).
#'
#' @inheritParams seq_refine
#' @inheritParams koutreads
#' @export
kractor_reads <- function(koutput, reads, ofile1 = NULL, ofile2 = NULL,
                          batch_size = NULL, chunk_bytes = NULL,
                          compression_level = 4L,
                          nqueue = NULL, threads = NULL, odir = NULL) {
    rust_kractor_reads(
        koutput = koutput,
        reads = reads,
        ofile1 = ofile1,
        ofile2 = ofile2,
        batch_size = batch_size,
        chunk_bytes = chunk_bytes,
        compression_level = compression_level,
        nqueue = nqueue,
        threads = threads,
        odir = odir
    )
}

rust_kractor_koutput <- function(kreport, koutput, ofile = NULL,
                                 taxonomy = c(
                                     "D__Bacteria", "D__Fungi", "D__Viruses"
                                 ),
                                 ranks = c("G", "S"),
                                 taxa = NULL,
                                 taxids = NULL,
                                 exclude = NULL,
                                 descendants = TRUE,
                                 batch_size = NULL, chunk_bytes = NULL,
                                 compression_level = 4L,
                                 nqueue = NULL, threads = NULL, odir = NULL,
                                 pprof = NULL) {
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    assert_string(ofile, allow_empty = FALSE, allow_null = TRUE)
    if (!is.null(taxonomy)) {
        taxonomy <- as.character(taxonomy)
        taxonomy <- taxonomy[!is.na(taxonomy)]
        if (length(taxonomy) == 0L) taxonomy <- NULL
    }
    if (!is.null(ranks)) {
        ranks <- as.character(ranks)
        ranks <- ranks[!is.na(ranks)]
        if (length(ranks) == 0L) ranks <- NULL
    }
    if (!is.null(taxa)) {
        taxa <- as.character(taxa)
        taxa <- taxa[!is.na(taxa)]
        if (length(taxa) == 0L) taxa <- NULL
    }
    if (!is.null(taxids)) {
        taxids <- as.character(taxids)
        taxids <- taxids[!is.na(taxids)]
        if (length(taxids) == 0L) taxids <- NULL
    }
    if (!is.null(exclude)) {
        exclude <- as.character(exclude)
        if (length(exclude) == 0L) exclude <- NULL
    }
    assert_bool(descendants)
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    assert_number_whole(chunk_bytes, min = 1, allow_null = TRUE)
    assert_number_whole(compression_level, min = 1, max = 12)
    assert_number_whole(threads,
        min = 0, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    threads <- threads %||% min(3, parallel::detectCores())
    nqueue <- check_queue(nqueue, 3L, threads)
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    assert_string(pprof, allow_empty = FALSE, allow_null = TRUE)
    odir <- odir %||% getwd()
    dir_create(odir)

    koutput_batch <- koutput_batch %||% KOUTPUT_BATCH
    fastq_batch <- fastq_batch %||% FASTQ_BATCH
    chunk_bytes <- chunk_bytes %||% CHUNK_BYTES

    ofile <- ofile %||% "kraken_microbiome_output.txt"
    ofile <- file.path(odir, ofile)

    if (is.null(pprof)) {
        rust_call(
            "kractor_koutput",
            kreport = kreport,
            koutput = koutput,
            taxonomy = taxonomy,
            ranks = ranks,
            taxa = taxa,
            taxids = taxids,
            exclude = exclude,
            descendants = descendants,
            ofile = ofile,
            compression_level = compression_level,
            batch_size = batch_size,
            chunk_bytes = chunk_bytes,
            nqueue = nqueue,
            threads = threads
        )
    } else {
        rust_call(
            "pprof_kractor_koutput",
            kreport = kreport,
            koutput = koutput,
            taxonomy = taxonomy,
            ranks = ranks,
            taxa = taxa,
            taxids = taxids,
            exclude = exclude,
            descendants = descendants,
            ofile = ofile,
            compression_level = compression_level,
            batch_size = batch_size,
            chunk_bytes = chunk_bytes,
            nqueue = nqueue,
            threads = threads,
            pprof_file = file.path(odir, pprof)
        )
    }
}

rust_kractor_reads <- function(koutput, reads, ofile1 = NULL, ofile2 = NULL,
                               batch_size = NULL, chunk_bytes = NULL,
                               compression_level = 4L,
                               nqueue = NULL, threads = NULL, odir = NULL,
                               pprof = NULL) {
    assert_string(koutput, allow_empty = FALSE)
    reads <- as.character(reads)
    if (length(reads) < 1L || length(reads) > 2L) {
        cli::cli_abort("{.arg reads} must be of length 1 or 2")
    }
    if (is_scalar(reads)) {
        fq1 <- reads[[1L]]
        fq2 <- NULL
    } else {
        fq1 <- reads[[1L]]
        fq2 <- reads[[2L]]
    }
    if ((is.null(fq2) && is.null(ofile1)) ||
        (!is.null(fq2) && is.null(ofile1) && is.null(ofile2))) {
        cli::cli_abort(c(
            "No output specified.",
            i = "Please provide at least one of {.arg ofile1} or {.arg ofile2} to write the results."
        ))
    }
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    assert_number_whole(chunk_bytes, min = 1, allow_null = TRUE)
    assert_number_whole(compression_level, min = 1, max = 12)
    assert_number_whole(threads,
        min = 0, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    threads <- threads %||% min(3, parallel::detectCores())
    nqueue <- check_queue(nqueue, 3L, threads)
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    assert_string(pprof, allow_empty = FALSE, allow_null = TRUE)
    odir <- odir %||% getwd()
    dir_create(odir)

    koutput_batch <- koutput_batch %||% KOUTPUT_BATCH
    fastq_batch <- fastq_batch %||% FASTQ_BATCH
    chunk_bytes <- chunk_bytes %||% CHUNK_BYTES

    if (is.null(pprof)) {
        rust_call(
            "kractor_reads",
            koutput = koutput,
            fq1 = fq1, ofile1 = file.path(odir, ofile1),
            fq2 = fq2, ofile2 = file.path(odir, ofile2),
            compression_level = compression_level,
            batch_size = batch_size,
            chunk_bytes = chunk_bytes,
            nqueue = nqueue,
            threads = threads
        )
    } else {
        rust_call(
            "pprof_kractor_reads",
            koutput = koutput,
            fq1 = fq1, ofile1 = file.path(odir, ofile1),
            fq2 = fq2, ofile2 = file.path(odir, ofile2),
            compression_level = compression_level,
            batch_size = batch_size,
            chunk_bytes = chunk_bytes,
            nqueue = nqueue,
            threads = threads,
            pprof_file = file.path(odir, pprof)
        )
    }
}

check_queue <- function(queue, default, threads, arg = caller_arg(queue),
                        call = rlang::caller_call()) {
    # styler: off
    if (.rlang_check_number(queue, min = 0,
                            allow_null = TRUE,
                            allow_decimal = FALSE,
                            allow_infinite = TRUE) != 0L ||
                            (!is.null(queue) && queue < 0)) {
        cli::cli_abort("{.arg {arg}} must be a non-negtive integer number")
    }
    # styler: on
    if (is.null(queue)) {
        default *
            if (is.null(threads) || threads == 0L) {
                parallel::detectCores()
            } else {
                threads
            }
    } else if (is.finite(queue)) {
        queue *
            if (is.null(threads) || threads == 0L) {
                parallel::detectCores()
            } else {
                threads
            }
    } else {
        NULL
    }
}
