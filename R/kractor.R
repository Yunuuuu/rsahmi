#' Extract Kraken2 Output and Reads (Kraken2 output and reads extractor)
#'
#' This function extracts reads and classification results from Kraken2 output,
#' based on a set of specified taxonomic IDs.
#'
#' @param kreport Path to the Kraken2 report file.
#' @param koutput Path to the Kraken2 output file.
#' @param reads A character vector of FASTQ files, either used as input to
#' Kraken2 or as classified output (recommended for efficiency due to smaller
#' size). Supports one file (single-end) or two files (paired-end).
#' @param ubread Path to the input sequence file that contains UMI and/or
#'   barcode sequences. If `NULL`, UMI/barcode parsing is disabled. When
#'   specified, `reads` must contain only a single FASTQ file. This is commonly
#'   used for 10x Genomics single-cell data, where the first read file contains
#'   only UMI and barcode sequences. In such cases, you should provide that file
#'   via `ubread`.
#' @param umi_ranges A range or a list of ranges specifying where UMI sequences
#'   are located in the reads. Must be created using the `ubrange()` function.
#'   Only used when `ubread` is not `NULL`.
#' @param barcode_ranges A range or a list of ranges specifying where cell
#'   barcode sequences are located in the reads. Must be created using the
#'   `ubrange()` function.  Only used when `ubread` is not `NULL`.
#' @param extract_koutput Path to the file where the extracted Kraken2 output
#' matching the specified `taxon` will be saved. Defaults to
#' `"kraken_microbiome_output.txt"`.
#' @param extract_reads A character vector of the same length as `reads`,
#' specifying the output FASTQ file(s) where the matching reads will be saved.
#' Defaults to `"kraken_microbiome_reads(_1|2).txt"`.
#' @param taxon An atomic character specify the taxa name wanted. Should follow
#' the kraken style, connected by rank codes, two underscores, and the
#' scientific name of the taxon (e.g., "d__Viruses").
#' @param batch_size Integer. Number of records to accumulate before triggering
#' a write operation. Default is `r code_quote(BATCH_SIZE, quote = FALSE)`.
#' @param chunk_size Integer. Size in bytes of the intermediate chunk used to
#'   split and distribute data to worker threads during processing. Default is
#'   `10 * 1024 * 1024` (10MB).
#' @param buffer_size Integer specifying the buffer size in bytes used for
#' writing to disk. This controls the capacity of the buffered file writer.
#' Default is `1 * 1024 * 1024` (1MB).
#' @param nqueue Integer. Maximum number of buffers per thread, controlling the
#'   amount of in-flight data awaiting writing. Default: `3`.
#' @param threads Integer. Number of threads to use. Default will determined
#' atomatically by rayon.
#' @param odir A string of directory to save the output files. Please see
#' `Value` section for details.
#' @param mmap Logical. Whether to enable memory-mapped file access. When set to
#'   `TRUE`, the function uses memory mapping, which can be highly efficient for
#'   multi-threaded reading and avoids redundant data copying. However, the
#'   performance of memory mapping may vary depending on the operating system
#'   and file system, and it is not always the fastest option. In most cases,
#'   standard file reading is already sufficiently fast. Therefore, the default
#'   is set to `FALSE`.
#' @seealso [`krakenx()`]
#' @return None. This function generates the following files:
#' - `extract_koutput`: Kraken2 output entries corresponding to the specified
#'   `taxon`, extracted from koutput.
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
                    ubread = NULL, umi_ranges = NULL, barcode_ranges = NULL,
                    extract_koutput = NULL, extract_reads = NULL,
                    taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                    chunk_size = NULL, buffer_size = NULL,
                    batch_size = NULL, nqueue = NULL,
                    threads = NULL, odir = NULL, mmap = FALSE) {
    rust_kractor_koutput(
        kreport, koutput,
        extract_koutput = extract_koutput,
        taxon = taxon,
        chunk_size = chunk_size,
        buffer_size = buffer_size,
        batch_size = batch_size,
        nqueue = nqueue,
        threads = threads,
        odir = odir,
        mmap = mmap
    )
    rust_kractor_reads(
        file.path(
            odir %||% getwd(),
            extract_koutput %||% "kraken_microbiome_output.txt"
        ),
        reads = reads,
        extract_reads = extract_reads,
        ubread = ubread,
        umi_ranges = umi_ranges,
        barcode_ranges = barcode_ranges,
        chunk_size = chunk_size,
        buffer_size = buffer_size,
        batch_size = batch_size,
        nqueue = nqueue,
        threads = threads,
        odir = odir,
        mmap = mmap
    )
    cli::cli_inform(c("v" = "Finished"))
}

BATCH_SIZE <- 500L
CHUNK_SIZE <- 10L * 1024L * 1024L
BUFFER_SIZE <- 1L * 1024L * 1024L

rust_kractor_koutput <- function(kreport, koutput, extract_koutput = NULL,
                                 taxon = c(
                                     "d__Bacteria", "d__Fungi",
                                     "d__Viruses"
                                 ),
                                 chunk_size = NULL, buffer_size = NULL,
                                 batch_size = NULL, nqueue = NULL,
                                 threads = NULL, odir = NULL,
                                 mmap = TRUE, pprof = NULL) {
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    assert_string(extract_koutput, allow_empty = FALSE, allow_null = TRUE)
    taxon <- as.character(taxon)
    if (length(taxon) == 0L) {
        cli::cli_abort("empty {.arg taxon} provided")
    }
    taxids <- parse_kraken_report(kreport)$
        explode(pl$col("ranks", "taxon"))$
        filter(
        pl$concat_str(
            pl$col("ranks")$str$to_lowercase(),
            pl$col("taxon"),
            separator = "__"
        )$is_in(pl$lit(taxon))
    )$
        select(pl$col("taxids")$list$last())$
        to_series()$unique()
    assert_number_whole(chunk_size, min = 1, allow_null = TRUE)
    assert_number_whole(buffer_size, min = 1, allow_null = TRUE)
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    assert_number_whole(threads,
        min = 0, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    if (is.null(odir)) odir <- getwd()
    assert_bool(mmap)
    assert_string(pprof, allow_empty = FALSE, allow_null = TRUE)
    dir_create(odir)

    # https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py#L95
    # take care of taxid: "A"
    if (taxids$is_in(pl$Series(values = c("81077", "A")))$any()) {
        taxids <- taxids$append(c("81077", "A"))
    }
    # the third column of kraken2 output:
    # Using (taxid ****)
    patterns <- paste0("(taxid ", as.character(taxids), ")")
    extract_koutput <- extract_koutput %||% "kraken_microbiome_output.txt"
    extract_koutput <- file.path(odir, extract_koutput)

    chunk_size <- chunk_size %||% CHUNK_SIZE
    buffer_size <- buffer_size %||% BUFFER_SIZE
    batch_size <- batch_size %||% BATCH_SIZE
    nqueue <- check_queue(nqueue, 3L) *
        if (is.null(threads) || threads == 0L) {
            parallel::detectCores()
        } else {
            threads
        }
    threads <- threads %||% 0L # Let rayon to determine the threads
    if (is.null(pprof)) {
        rust_call(
            "kractor_koutput",
            patterns = patterns,
            koutput = koutput,
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
            patterns = patterns,
            koutput = koutput,
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

rust_kractor_reads <- function(koutput, reads,
                               ubread = NULL,
                               umi_ranges = NULL, barcode_ranges = NULL,
                               extract_reads = NULL,
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
    assert_string(ubread, allow_empty = FALSE, allow_null = TRUE)
    if (!is.null(umi_ranges) && !is_range(umi_ranges)) {
        cli::cli_abort("{.arg umi_ranges} must be created with {.fn ubrange} or a combination of them using {.fn c}.")
    }
    if (!is.null(barcode_ranges) && !is_range(barcode_ranges)) {
        cli::cli_abort("{.arg barcode_ranges} must be created with {.fn ubrange} or a combination of them using {.fn c}.")
    }
    if (!identical(is.null(ubread), is.null(umi_ranges)) ||
        !identical(is.null(ubread), is.null(barcode_ranges))) {
        cli::cli_abort(c(
            "All or none of the following arguments must be provided:",
            "*" = "{.arg ubread}",
            "*" = "{.arg umi_ranges}",
            "*" = "{.arg barcode_ranges}",
            "i" = "These are required to extract {.field UMI} and {.field Cell Barcode}."
        ))
    }
    if (!is.null(ubread) && length(reads) == 2L) {
        cli::cli_abort(c(
            "{.arg reads} must be a single FASTQ file when {.arg ubread} is specified.",
            i = paste(
                "UMI/barcode parsing requires {.arg reads} to contain only one input file",
                "while the read containing UMI/barcode sequences must be passed via {.arg ubread}."
            )
        ))
    }
    assert_number_whole(chunk_size, min = 1, allow_null = TRUE)
    assert_number_whole(buffer_size, min = 1, allow_null = TRUE)
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    assert_number_whole(threads,
        min = 1, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    if (is.null(odir)) odir <- getwd()
    assert_bool(mmap)
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
    nqueue <- check_queue(nqueue, 3L) *
        if (is.null(threads) || threads == 0L) {
            parallel::detectCores()
        } else {
            threads
        }
    threads <- threads %||% 0L
    if (is.null(pprof)) {
        rust_call(
            "kractor_reads",
            koutput,
            fq1 = fq1, ofile1 = extract_read1,
            fq2 = fq2, ofile2 = extract_read2,
            ubread = ubread, umi_ranges = umi_ranges,
            barcode_ranges = barcode_ranges,
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
            ubread = ubread, umi_ranges = umi_ranges,
            barcode_ranges = barcode_ranges,
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

check_queue <- function(queue, default, arg = caller_arg(queue),
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
        default
    } else if (is.finite(queue)) {
        queue
    } else {
        NULL
    }
}
