#' Extract reads and output from Kraken
#'
#' @param kraken_report The path to kraken report file.
#' @param reads The original fastq files (input in `kraken2`). You can pass
#' two paired-end files directly.
#' @param ...
#'  - `extract_kraken_output`: Additional arguments passed to
#'    [sink_csv][polars::LazyFrame_sink_csv].
#'  - `extract_kraken_reads`: Additional arguments passed to
#'    [run][biosys::Execute] method.
#' @name extractor
#' @seealso
#' <https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown>
NULL

#' @param taxon An atomic character specify the taxa name wanted. Should follow
#' the kraken style, connected by rank codes, two underscores, and the
#' scientific name of the taxon (e.g., "d__Viruses")
#' @export
#' @rdname extractor
#' @importFrom polars pl
extract_taxids <- function(kraken_report, taxon = c("d__Bacteria", "d__Fungi", "d__Viruses")) {
    parse_kraken_report(kraken_report)$
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
}

#' @param kraken_out The path to kraken output file.
#' @param taxids A character specify NCBI taxonony identifier to extract.
#' @param ofile A string of file save the kraken output of specified `taxids`.
#' @param odir A string of directory to save the `ofile`.
#' @export
#' @rdname extractor
#' @importFrom polars pl
extract_kraken_output <- function(kraken_out, taxids,
                                  ofile = "kraken_microbiome_output.txt",
                                  odir = getwd(), ...) {
    dir_create(odir)
    # https://github.com/jenniferlu717/KrakenTools/blob/master/extract_kraken_reads.py#L95
    # take care of taxid: "A"
    if (taxids$is_in(pl$Series(values = c("81077", "A")))$any()) {
        taxids <- taxids$append(c("81077", "A"))
    }
    pl$scan_csv(kraken_out, has_header = FALSE, separator = "\t")$
        filter(
        pl$col("column_3")$str$
            extract(pl$lit("\\s*(.+)\\s*\\(taxid\\s*(\\d+|A)\\s*\\)"), 2L)$
            is_in(taxids)
    )$
        sink_csv(
        path = file.path(odir, ofile),
        include_header = FALSE,
        separator = "\t", ...
    )
}

#' @param threads Number of threads to use, see `biosys::seqkit("grep")$help()`.
#' @inheritParams biosys::seqkit
#' @export
#' @rdname extractor
#' @importFrom polars pl
extract_kraken_reads <- function(kraken_out, reads, ofile = NULL,
                                 odir = getwd(), threads = NULL,
                                 ..., seqkit = NULL) {
    if (length(reads) < 1L || length(reads) > 2L) {
        cli::cli_abort("{.arg reads} must be of length 1 or 2")
    }
    if (is.null(ofile)) {
        if (is_scalar(reads)) {
            ofile <- "kraken_microbiome_reads.fa"
        } else {
            ofile <- sprintf("kraken_microbiome_reads_%d.fa", seq_along(reads))
        }
    } else if (length(ofile) != length(reads)) {
        cli::cli_abort("{.arg ofile} must have the same length of {.arg reads}")
    }
    dir_create(odir)
    ofile <- file.path(odir, ofile)
    file <- tempfile("kraken_sequence_id")
    pl$scan_csv(kraken_out, has_header = FALSE, separator = "\t")$
        # second column is the sequence id
        select(pl$col("column_2"))$unique()$
        sink_csv(path = file, include_header = FALSE, separator = "\t")
    on.exit(file.remove(file))
    status <- vapply(seq_along(reads), function(i) {
        extract_sequence_id(
            fq = reads[[i]], ofile = ofile[[i]],
            sequence_id = file, ...,
            threads = threads, seqkit = seqkit
        )
    }, integer(1L))
    invisible(status)
}

extract_sequence_id <- function(fq, ofile, sequence_id, ..., threads, seqkit) {
    biosys::seqkit("seq", "--only-id", fq, seqkit = seqkit)$
        pipe(
        biosys::seqkit(
            "grep", biosys::arg("-f", sequence_id),
            if (!is.null(threads)) {
                biosys::arg("--threads", threads, format = "%d")
            },
            "-n",
            seqkit = seqkit
        )
    )$
        pipe(
        biosys::seqkit("fq2fa", biosys::arg("-o", ofile), seqkit = seqkit)
    )$
        run(...)
}
