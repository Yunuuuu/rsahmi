#' Extract Kraken2 Output and Reads (Kraken2 output and reads extractor)
#'
#' This function extracts reads and classification results from Kraken2 output,
#' based on a set of specified taxonomic IDs.
#'
#' @param kreport Path to the Kraken2 report file.
#' @param koutput Path to the Kraken2 output file.
#' @param reads A character vector of FASTQ file paths, either the original
#' reads used as input to Kraken2 or the classified output reads (recommended
#' for efficiency as they are smaller). Accepts one file for single-end or two
#' files for paired-end. Note: if `ubread` is specified, `reads` must be the
#' original Kraken2 input reads to ensure correct barcode/UMI pairing.
#' @param extract_koutput Path to the file where the extracted Kraken2 output
#' matching the specified `taxon` will be saved. Defaults to
#' `"kraken_microbiome_output.txt"`.
#' @param extract_reads A character vector of the same length as `reads`,
#' specifying the output FASTQ file(s) where the matching reads will be saved.
#' Defaults to `"kraken_microbiome_reads(_1|2).txt"`.
#' @param taxonomy Character vector. The set of taxonomic groups to include
#' (default: `c("d__Bacteria", "d__Fungi", "d__Viruses")`). This defines the
#' global taxa to consider. If `NULL`, all taxa will be used. If `descendants =
#' TRUE`, only the descendants within these groups will be considered. The
#' selection of taxa can be further refined using the `ranks`, `taxa`, and
#' `taxids` parameters. One of `taxonomy`, `ranks`, `taxa`, or `taxids` must be
#' provided.
#' @param ranks Character vector. The taxonomic ranks to filter by. Default:
#' `c("G", "S")` (optional).
#' @param taxa Character vector. Specific taxa to include (optional).
#' @param taxids Character vector. A list of taxid values to filter by
#' (optional).
#' @param descendants Logical. Whether to include descendants of the selected
#' taxa (default: `TRUE`).
#' @inheritParams seq_refine
#' @param mmap_koutput,mmap_reads Logical. Whether to enable memory-mapped file
#'   access for reading input FASTQ reads (`mmap_reads`) or Kraken2 output
#'   (`mmap_koutput`). When set to `TRUE`, the function uses memory mapping to
#'   reduce data copying and improve performance in multi-threaded environments.
#'   This can be highly efficient on some systems, but performance gains may
#'   vary depending on the operating system and file system. Defaults to `TRUE`.
#' @seealso [`krakenx()`]
#' @return None. This function generates the following files:
#' - `extract_koutput`: Kraken2 output entries corresponding to the specified
#'   `taxonomy`, `ranks`, `taxa`, `taxids`, and `descendants` extracted from the
#'   input `koutput`.
#' - `extract_reads`: Sequence file(s) containing reads assigned to the
#'   specified `taxon`.
#' @examples
#' \dontrun{
#' blit::kraken2(
#'     fq1 = fq1,
#'     fq2 = fq2,
#'     classified_out = "classified.fq",
#'     # Number of threads to use
#'     blit::arg("--threads", 10L, format = "%d"),
#'     # the kraken database
#'     blit::arg("--db", kraken_db),
#'     "--use-names", "--report-minimizer-data",
#' ) |> blit::cmd_run()
#'
#' # 1. `kreport` should be the kraken2 report file of `blit::kraken2()`
#' # 2. `koutput` should be the kraken2 output file of `blit::kraken2()`
#' # 3. `reads` should be the same with `fq1` and `fq2` in `blit::kraken2()`
#' kractor(
#'     kreport = "kraken_report.txt",
#'     koutput = "kraken_output.txt",
#'     reads = c(fq1, fq2)
#' )
#' }
#' @export
kractor <- function(kreport, koutput, reads,
                    extract_koutput = NULL, extract_reads = NULL,
                    taxonomy = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                    ranks = c("G", "S"),
                    taxa = NULL,
                    taxids = NULL,
                    descendants = TRUE,
                    chunk_size = NULL, buffer_size = NULL,
                    batch_size = NULL, nqueue = NULL,
                    threads = NULL, odir = NULL,
                    mmap_koutput = TRUE, mmap_reads = TRUE) {
    rust_kractor_koutput(
        kreport = kreport,
        koutput = koutput,
        extract_koutput = extract_koutput,
        taxonomy = taxonomy,
        ranks = ranks,
        taxa = taxa,
        taxids = taxids,
        descendants = descendants,
        chunk_size = chunk_size,
        buffer_size = buffer_size,
        batch_size = batch_size,
        nqueue = nqueue,
        threads = threads,
        odir = odir,
        mmap = mmap_koutput
    )
    rust_kractor_reads(
        koutput = file.path(
            odir %||% getwd(),
            extract_koutput %||% "kraken_microbiome_output.txt"
        ),
        reads = reads,
        extract_reads = extract_reads,
        chunk_size = chunk_size,
        buffer_size = buffer_size,
        batch_size = batch_size,
        nqueue = nqueue,
        threads = threads,
        odir = odir,
        mmap = mmap_reads
    )
    cli::cli_inform(c("v" = "Finished"))
}

BATCH_SIZE <- 500L
FASTQ_BATCH <- 256
KOUTPUT_BATCH <- 1000
CHUNK_BYTES <- 10L * 1024L * 1024L

#' Filter Kraken2 Output by Taxon
#'
#' This function filters Kraken2 classification output (`koutput`) by taxonomic
#' identifiers derived from a Kraken2 report (`kreport`). It extracts lines
#' matching the desired `taxonomy`, `ranks`, `taxa`, `taxids`, and `descendants`
#' and writes the filtered results to an output file.
#'
#' @inheritParams kractor
#' @param mmap Logical. Whether to enable memory-mapped file access for reading
#'   input Kraken2 output. When set to `TRUE`, the function uses memory mapping
#'   to reduce data copying and improve performance in multi-threaded
#'   environments. This can be highly efficient on some systems, but
#'   performance gains may vary depending on the operating system and file
#'   system. Defaults to `TRUE`.
#' @return None. The function generates a filtered Kraken2 output file
#'   containing entries corresponding to the specified `taxonomy`, `ranks`,
#'   `taxa`, `taxids`, and `descendants` extracted from the input `koutput`.
#' @export
kractor_koutput <- function(kreport, koutput, extract_koutput = NULL,
                            taxonomy = c(
                                "d__Bacteria", "d__Fungi",
                                "d__Viruses"
                            ),
                            ranks = c("G", "S"),
                            taxa = NULL,
                            taxids = NULL,
                            descendants = TRUE,
                            chunk_size = NULL, buffer_size = NULL,
                            batch_size = NULL, nqueue = NULL,
                            mmap = TRUE, threads = NULL, odir = NULL) {
    rust_kractor_koutput(
        kreport = kreport,
        koutput = koutput,
        extract_koutput = extract_koutput,
        taxonomy = taxonomy,
        ranks = ranks,
        taxa = taxa,
        taxids = taxids,
        descendants = descendants,
        chunk_size = chunk_size,
        buffer_size = buffer_size,
        batch_size = batch_size,
        nqueue = nqueue,
        threads = threads,
        odir = odir,
        mmap = mmap
    )
}

#' Extract Reads from Kraken2 Output Based on Classification
#'
#' This function extracts reads corresponding to selected classifications from a
#' Kraken2 output file (`koutput`). Only reads classified to selected taxa will
#' be extracted from the provided sequence file (`reads`).
#'
#' @inheritParams kractor
#' @param mmap Logical. Whether to enable memory-mapped file access for reading
#'   input sequence files. When set to `TRUE`, the function uses memory mapping
#'   to reduce data copying and improve performance in multi-threaded
#'   environments. This can be highly efficient on some systems, but
#'   performance gains may vary depending on the operating system and file
#'   system. Defaults to `TRUE`.
#' @export
kractor_reads <- function(koutput, reads, extract_reads = NULL,
                          chunk_size = NULL, buffer_size = NULL,
                          batch_size = NULL,
                          nqueue = NULL, threads = NULL,
                          odir = NULL, mmap = TRUE) {
    rust_kractor_reads(
        koutput = koutput,
        reads = reads,
        extract_reads = extract_reads,
        chunk_size = chunk_size,
        buffer_size = buffer_size,
        batch_size = batch_size,
        nqueue = nqueue,
        threads = threads,
        odir = odir,
        mmap = mmap
    )
}

rust_kractor_koutput <- function(kreport, koutput, extract_koutput = NULL,
                                 taxonomy = c(
                                     "d__Bacteria", "d__Fungi",
                                     "d__Viruses"
                                 ),
                                 ranks = c("G", "S"),
                                 taxa = NULL,
                                 taxids = NULL,
                                 descendants = TRUE,
                                 chunk_size = NULL, buffer_size = NULL,
                                 batch_size = NULL, nqueue = NULL,
                                 threads = NULL, odir = NULL, mmap = TRUE,
                                 pprof = NULL) {
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    assert_string(extract_koutput, allow_empty = FALSE, allow_null = TRUE)
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
    assert_bool(descendants)
    assert_number_whole(chunk_size, min = 1, allow_null = TRUE)
    assert_number_whole(buffer_size, min = 1, allow_null = TRUE)
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    assert_number_whole(threads,
        min = 0, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    nqueue <- check_queue(nqueue, 3L, threads)
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    assert_bool(mmap)
    if (is.null(odir)) odir <- getwd()
    assert_string(pprof, allow_empty = FALSE, allow_null = TRUE)
    dir_create(odir)

    extract_koutput <- extract_koutput %||% "kraken_microbiome_output.txt"
    extract_koutput <- file.path(odir, extract_koutput)

    chunk_size <- chunk_size %||% CHUNK_SIZE
    buffer_size <- buffer_size %||% BUFFER_SIZE
    batch_size <- batch_size %||% BATCH_SIZE
    threads <- threads %||% 0L # Let rayon to determine the threads
    if (is.null(pprof)) {
        rust_call(
            "kractor_koutput",
            kreport = kreport,
            koutput = koutput,
            taxonomy = taxonomy,
            ranks = ranks,
            taxa = taxa,
            taxids = taxids,
            descendants = descendants,
            ofile = extract_koutput,
            chunk_size = chunk_size,
            buffer_size = buffer_size,
            batch_size = batch_size,
            nqueue = nqueue,
            threads = threads,
            mmap = mmap
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
            descendants = descendants,
            ofile = extract_koutput,
            chunk_size = chunk_size,
            buffer_size = buffer_size,
            batch_size = batch_size,
            nqueue = nqueue,
            threads = threads,
            mmap = mmap,
            pprof_file = file.path(odir, pprof)
        )
    }
}

rust_kractor_reads <- function(koutput, reads, extract_reads = NULL,
                               chunk_size = NULL, buffer_size = NULL,
                               batch_size = NULL,
                               nqueue = NULL, threads = NULL,
                               odir = NULL, mmap = TRUE, pprof = NULL) {
    assert_string(koutput, allow_empty = FALSE)
    reads <- as.character(reads)
    if (length(reads) < 1L || length(reads) > 2L) {
        cli::cli_abort("{.arg reads} must be of length 1 or 2")
    }
    if (is.null(extract_reads)) {
        if (is_scalar(reads)) {
            extract_reads <- "kraken_microbiome_reads.fa"
        } else {
            extract_reads <- sprintf(
                "kraken_microbiome_reads_%d.fa",
                seq_along(reads)
            )
        }
    } else {
        extract_reads <- as.character(extract_reads)
        if (length(extract_reads) != length(reads)) {
            cli::cli_abort(
                "{.arg extract_reads} must have the same length of {.arg reads}"
            )
        }
    }
    assert_number_whole(chunk_size, min = 1, allow_null = TRUE)
    assert_number_whole(buffer_size, min = 1, allow_null = TRUE)
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    assert_number_whole(threads,
        min = 1, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    nqueue <- check_queue(nqueue, 3L, threads)
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    assert_bool(mmap)
    if (is.null(odir)) odir <- getwd()
    assert_string(pprof, allow_empty = FALSE, allow_null = TRUE)
    dir_create(odir)

    extract_reads <- file.path(odir, extract_reads)
    if (is_scalar(extract_reads)) {
        fq1 <- reads[[1L]]
        fq2 <- NULL
        extract_read1 <- extract_reads[[1L]]
        extract_read2 <- NULL
    } else {
        fq1 <- reads[[1L]]
        fq2 <- reads[[2L]]
        extract_read1 <- extract_reads[[1L]]
        extract_read2 <- extract_reads[[2L]]
    }

    chunk_size <- chunk_size %||% CHUNK_SIZE
    buffer_size <- buffer_size %||% BUFFER_SIZE
    batch_size <- batch_size %||% BATCH_SIZE
    threads <- threads %||% 0L
    if (is.null(pprof)) {
        rust_call(
            "kractor_reads",
            koutput,
            fq1 = fq1, ofile1 = extract_read1,
            fq2 = fq2, ofile2 = extract_read2,
            chunk_size = chunk_size,
            buffer_size = buffer_size,
            batch_size = batch_size,
            nqueue = nqueue,
            threads = threads,
            mmap = mmap
        )
    } else {
        rust_call(
            "pprof_kractor_reads",
            koutput,
            fq1 = fq1, ofile1 = extract_read1,
            fq2 = fq2, ofile2 = extract_read2,
            chunk_size = chunk_size,
            buffer_size = buffer_size,
            batch_size = batch_size,
            nqueue = nqueue,
            threads = threads,
            mmap = mmap,
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
