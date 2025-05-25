#' Extract Reads and Output from Kraken2 (kraken2 extractor)
#'
#' This function extracts reads and classification results from Kraken2 output,
#' based on a set of specified taxonomic IDs.
#'
#' @param kreport Path to the Kraken2 report file.
#'
#' @param koutput Path to the Kraken2 output file.
#'
#' @param reads A character vector of FASTQ files used as input to Kraken2.
#' Can be one file (single-end) or two files (paired-end).
#'
#' @param extract_koutput Path to the file where the extracted Kraken2 output
#' matching the specified `taxon` will be saved. Defaults to
#' `"kraken_microbiome_output.txt"`.
#'
#' @param extract_reads A character vector of the same length as `reads`,
#' specifying the output FASTQ file(s) where the matching reads will be saved.
#' Defaults to `"kraken_microbiome_reads(_1|2).txt"`
#'
#' @param taxon An atomic character specify the taxa name wanted. Should follow
#' the kraken style, connected by rank codes, two underscores, and the
#' scientific name of the taxon (e.g., "d__Viruses").
#'
#' @param write_buffer Integer specifying the buffer size in bytes used for
#' writing to disk. This controls the capacity of the buffered file writer.
#' Default is `1 * 1024 * 1024` (1MB).
#'
#' @param read_buffer Integer specifying the size in bytes of the intermediate
#' buffer used for splitting and distributing chunks to worker threads during
#' processing. Default is `2 * 1024 * 1024` (2MB).
#'
#' @param parse_buffer Integer. Number of records to write per batch. Default is
#'   `10`.
#'
#' @param read_queue,write_queue Integer. Maximum number of buffers per thread,
#'   controlling the amount of in-flight data awaiting processing or
#'   writing. Default is `2`.
#'
#' @param threads Integer. Number of threads to use. Defaults to all available
#' threads.
#'
#' @param odir A string of directory to save the `ofile`.
#'
#' @seealso [`kraken_taxon()`]
#' @return None. This function generates the following files:
#' - `extract_koutput`: Kraken2 output entries corresponding to the specified
#'   `taxon`, extracted from koutput.
#' - `extract_reads`: Sequence file(s) containing reads assigned to the
#'   specified `taxon`.
#' @examples
#' \dontrun{
#' # For 10x Genomic data, `fq1` only contain barcode and umi, but the official
#' # didn't give any information for this. In this way, I prefer using
#' # `umi-tools` to transform the `umi` into fq2 and then run `rsahmi` with
#' # only fq2.
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
                    taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                    read_buffer = NULL, write_buffer = NULL,
                    parse_buffer = NULL,
                    read_queue = NULL, write_queue = NULL,
                    threads = NULL, odir = getwd()) {
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    reads <- as.character(reads)
    if (length(reads) < 1L || length(reads) > 2L) {
        cli::cli_abort("{.arg reads} must be of length 1 or 2")
    }
    assert_string(extract_koutput, allow_empty = FALSE, allow_null = TRUE)
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
    extract_koutput <- extract_koutput %||% "kraken_microbiome_output.txt"
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
    assert_string(odir, allow_empty = FALSE)
    dir_create(odir)

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

    rust_kractor(
        koutput, reads, taxids,
        read_buffer, write_buffer, parse_buffer,
        read_queue, write_queue,
        extract_reads, extract_koutput, threads, odir
    )
    cli::cli_inform(c("v" = "Finished"))
}

rust_kractor <- function(koutput, reads, taxids,
                         read_buffer, write_buffer, parse_buffer,
                         read_queue, write_queue,
                         extract_reads, extract_koutput,
                         threads, odir) {
    # https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py#L95
    # take care of taxid: "A"
    if (taxids$is_in(pl$Series(values = c("81077", "A")))$any()) {
        taxids <- taxids$append(c("81077", "A"))
    }
    # the third column of kraken2 output:
    # Using (taxid ****)
    patterns <- paste0("(taxid ", as.character(taxids), ")")
    read_buffer <- read_buffer %||% (2 * 1024L * 1024L) # DEFAULT_BUF_SIZE 2MB
    write_buffer <- write_buffer %||% (1 * 1024L * 1024L) # 1MB
    parse_buffer <- parse_buffer %||% 10L
    read_queue <- read_queue %||% 2L
    write_queue <- write_queue %||% 2L
    extract_koutput <- file.path(odir, extract_koutput)
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
    rust_call(
        "kractor", koutput, patterns,
        ofile = extract_koutput,
        fq1 = fq1, ofile1 = extract_read1,
        fq2 = fq2, ofile2 = extract_read2,
        read_buffer = read_buffer,
        write_buffer = write_buffer,
        parse_buffer = parse_buffer,
        read_queue = read_queue,
        write_queue = write_queue,
        threads = threads
    )
}
