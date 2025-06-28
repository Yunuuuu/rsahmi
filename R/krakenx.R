#' Full SAHMI-Compatible Pipeline: Kraken2 + Read Extraction + K-mer/UMI
#' Quantification
#'
#' @inheritParams blit::kraken2
#' @param db Path to the Kraken2 database. You can download prebuilt databases
#'   from <https://benlangmead.github.io/aws-indexes/k2>, or build your own by
#'   following the instructions at
#'   <https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases>.
#' @inheritParams kractor
#' @inheritParams kuactor
#' @param kraken2_threads Integer. Number of threads to use. Default to all
#' available threads.
#' @param kractor_threads Integer. Number of threads to use. Default will
#' determined atomatically by rayon.
#' @param kuactor_threads Integer. Number of threads to use. Default to all
#' available threads.
#' @param envpath String. Additional path to include in the `PATH` environment
#'   variable. Used to locate the `kraken2` executable.
#' @param overwrite Logical. Whether to overwrite existing files in `odir`.
#' @seealso
#' - [`kractor()`]
#' - [`kuactor()`]
#' @return None. This function generates the following files:
#' - `kreport`: the `kraken2` report file.
#' - `koutput`: the `kraken2` output file.
#' - `classified_out`: the classified sequence file(s) from `kraken2`. Not
#'   required for downstream processing.
#' - `unclassified_out`: the unclassified sequence file(s) from `kraken2`. Not
#'   required for downstream processing.
#' - `extract_koutput`: Kraken2 output entries corresponding to the specified
#'   `taxon`, extracted from `koutput`.
#' - `extract_reads`: Sequence file(s) containing reads assigned to the
#'   specified `taxon`.
#'  - `extract_kreport`: A filtered version of the Kraken2 taxonomic report,
#'   containing only taxa that meet the `ranks` criteria and are observed in
#'   `koutput`. Used by [`rpmm_quantile()`] and [`slsd()`].
#'  - `extract_kmer`: A table quantifying total and unique k-mers assigned to
#'   each taxon across barcodes. Used by [`blsd()`].
#'  - `extract_umi`: A table of taxon–barcode–UMI combinations indicating all
#'   observed UMI-tagged reads per taxon. Used by [`taxa_counts()`].
#' @export
krakenx <- function(reads, ...,
                    ubread = NULL, umi_ranges = NULL, barcode_ranges = NULL,
                    db = NULL,
                    kreport = "kraken_report.txt",
                    koutput = "kraken_output.txt",
                    classified_out = "classified.fq",
                    unclassified_out = NULL,
                    extract_koutput = NULL,
                    extract_reads = NULL,
                    extract_kreport = NULL,
                    extract_kmer = NULL,
                    extract_umi = NULL,
                    taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                    ranks = c("G", "S"),
                    chunk_size = NULL, buffer_size = NULL,
                    batch_size = NULL, nqueue = NULL,
                    kmer_len = 35L, min_frac = 0.5, exclude = "9606",
                    barcode_extractor = function(sequence_id, read1, read2) {
                        substring(read1, 1L, 16L)
                    },
                    umi_extractor = function(sequence_id, read1, read2) {
                        substring(read1, 17L, 28L)
                    },
                    mmap = FALSE,
                    kraken2_threads = NULL,
                    kractor_threads = NULL,
                    kuactor_threads = NULL,
                    kraken2 = NULL, envpath = NULL,
                    overwrite = FALSE, odir = getwd()) {
    assert_bool(overwrite)
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    if (is.null(odir)) odir <- getwd()

    # For kraken2 command
    reads <- as.character(reads)
    if (length(reads) == 1L) {
        fq1 <- reads
        fq2 <- NULL
    } else if (length(reads) > 2L) {
        fq1 <- reads[1L]
        fq2 <- reads[2L]
    } else {
        cli::cli_abort("{.arg reads} must be of length 1 or 2")
    }
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    assert_string(classified_out, allow_empty = FALSE)
    assert_string(unclassified_out, allow_empty = FALSE, allow_null = TRUE)
    kraken2_files <- file.path(
        odir,
        c(koutput, kreport, classified_out, unclassified_out)
    )
    if (overwrite || !all(file.exists(kraken2_files))) {
        assert_string(db, allow_empty = FALSE, allow_null = TRUE)
        assert_number_whole(kraken2_threads,
            min = 1, max = as.double(parallel::detectCores()),
            allow_null = TRUE
        )
        if (!is.null(db)) db <- sprintf("--db %s", db)
        kraken2_threads <- kraken2_threads %||% parallel::detectCores()
        command <- blit::kraken2(
            fq1 = fq1,
            fq2 = fq2,
            ...,
            ofile = koutput,
            report = kreport,
            classified_out = classified_out,
            unclassified_out = unclassified_out,
            db,
            "--use-names",
            "--report-minimizer-data",
            sprintf("--threads %d", kraken2_threads),
            odir = odir,
            kraken2 = kraken2
        )
        command <- blit::cmd_envpath(command, envpath)
        blit::cmd_run(command, spinner = TRUE, verbose = TRUE)
    } else {
        cli::cli_inform("Using files exist: {.path {kraken2_files}}")
    }
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
                "sequence file {.arg reads}"
            ))
        }
    }
    extract_koutput <- extract_koutput %||% "kraken_microbiome_output.txt"
    kractor_files <- file.path(odir, c(extract_koutput, extract_reads))
    if (overwrite || !all(file.exists(kractor_files))) {
        assert_number_whole(kractor_threads,
            min = 1, max = as.double(parallel::detectCores()),
            allow_null = TRUE
        )
        kractor(
            kreport = file.path(odir, kreport),
            koutput = file.path(odir, koutput),
            reads = reads, ubread = ubread,
            umi_ranges = umi_ranges, barcode_ranges = barcode_ranges,
            extract_koutput = extract_koutput,
            extract_reads = extract_reads,
            taxon = taxon,
            chunk_size = chunk_size,
            buffer_size = buffer_size,
            batch_size = batch_size,
            nqueue = nqueue,
            threads = kractor_threads,
            odir = odir,
            mmap = mmap
        )
    } else {
        cli::cli_inform("Using file exist: {.path {kractor_files}}")
    }

    extract_kreport <- extract_kreport %||% "kraken_microbiome_report.txt"
    extract_kmer <- extract_kmer %||% "kraken_microbiome_kmer.txt"
    extract_umi <- extract_umi %||% "kraken_microbiome_umi.txt"
    kuactor_files <- file.path(
        odir,
        c(extract_kreport, extract_kmer, extract_umi)
    )
    if (overwrite || !all(file.exists(kuactor_files))) {
        assert_number_whole(kuactor_threads,
            min = 1, max = as.double(parallel::detectCores()),
            allow_null = TRUE
        )
        kuactor(
            kreport = file.path(odir, kreport),
            koutput = file.path(odir, extract_koutput),
            reads = file.path(odir, extract_reads),
            extract_kreport = extract_kreport,
            extract_kmer = extract_kmer,
            extract_umi = extract_umi,
            barcode_extractor = barcode_extractor,
            umi_extractor = umi_extractor,
            ranks = ranks,
            kmer_len = kmer_len,
            min_frac = min_frac,
            exclude = exclude,
            threads = kuactor_threads,
            odir = odir
        )
    } else {
        cli::cli_inform("Using file exist: {.path {kuactor_files}}")
    }

    if (!overwrite &&
        all(file.exists(c(kraken2_files, kractor_files, kuactor_files)))) {
        cli::cli_warn(paste(
            "Nothing to do, please set {.code overwrite = TRUE}",
            "if you want to overwrite the files exist"
        ))
    }
    invisible(NULL)
}
