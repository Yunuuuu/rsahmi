#' Quantitation of microbes
#' @description After identifying true taxa, reads assigned to those taxa are
#' extracted and then undergo a series of filters. The cell barcode and UMI are
#' used to demultiplex the reads and create a barcode x taxa counts matrix. The
#' full taxonomic classification of all resulting barcodes and the number of
#' counts assigned to each clade are tabulated.
#' @param taxa Taxa ids to extracted, can be a file, in which, the first column
#'   will be extracted.
#' @param n_filter Filter reads with `> n_filter` of one nucleotide.
#' @inheritParams run_sckmer
#' @importFrom parallelly availableCores
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @export
taxa_counts <- function(fa1, fa2, kraken_report, mpa_report, taxa, sample = NULL, out_dir = getwd(), cb_len = 16L, umi_len = 10L, n_filter = 130L, cores = availableCores()) {
    # read in Fasta data -----------------------------------------
    reads <- ShortRead::readFasta(fa1)

    # read in taxa data ------------------------------------------
    if (is.character(taxa) && length(taxa) == 1L && file.exists(taxa)) {
        taxa <- data.table::fread(taxa, header = FALSE)[[1L]]
    }
    taxa <- as.integer(taxa)
    taxa <- taxa[!is.na(taxa)]
    # taxa_list <- split(taxa, ceiling(seq_along(taxa) / 2000L))
    # taxa_list <- lapply(taxa_list, function(tx) {
    #     paste0("\\*", tx, "\\*", collapse = "|")
    # })

    # filter reads -----------------------------------------------
    filter_reads <- ShortRead::polynFilter(
        threshold = n_filter, nuc = c("A", "T", "G", "C")
    )
    filter_reads <- ShortRead::compose(filter_reads)
    reads <- reads[filter_reads(reads)]

    # extract necessary infos ------------------------------------
    sequences <- ShortRead::sread(reads)
    headers <- ShortRead::id(reads)
    barcode <- substr(sequences, 1L, cb_len)
    umi <- substr(sequences, cb_len + 1L, cb_len + umi_len)
    taxid <- gsub(".*taxid\\|", "", headers, perl = TRUE)
    length(barcode)

    # prepare kraken report and mpa report data ------------------
    kr <- data.table::fread(kraken_report, header = FALSE, sep = "\t")[-c(1:2)]
    kr[, V8 := gsub("[^[:alnum:]]+", "_", V8, perl = TRUE)] # nolint

    mpa <- data.table::fread(mpa_report, header = FALSE, sep = "\t")
    mpa[, taxid := vapply(strsplit(V1, "|", fixed = TRUE), function(x) { # nolint
        str <- sub(
            ".*__", "",
            unlist(x, recursive = FALSE, use.names = FALSE),
            perl = TRUE
        )
        str <- gsub("[^[:alnum:]]+", "_", str, perl = TRUE)
        paste0("*", paste0(
            kr$V7[data.table::chmatch(str, kr$V8)],
            collapse = "*"
        ), "*")
    }, character(1L))]
    # mpa$taxid[1L] <- NA_character_

    # extract taxa data -----------------------------------------
    full_taxa <- lapply(strsplit(mpa$taxid, "*", fixed = TRUE), function(x) {
        x <- x[x != "NA" & x != "" & !is.na(x)]
        as.integer(x)
    })
    cli::cli_alert("Finding children taxa")
    old_handlers <- new_handlers()
    on.exit(progressr::handlers(old_handlers))
    old_plan <- future::plan("multicore", workers = cores)
    on.exit(future::plan(old_plan), add = TRUE)
    p <- progressr::progressor(along = full_taxa, auto_finish = TRUE)
    child_taxa <- future.apply::future_lapply(full_taxa, function(x) {
        p()
        x[cumsum(x %in% taxa) > 0L]
    })
    child_taxa <- unlist(child_taxa, recursive = FALSE, use.names = FALSE)

    # create barcode umi data -----------------------------------
    idx <- which(taxid %in% child_taxa)
    barcode <- barcode[idx]
    umi <- umi[idx]
    taxid <- taxid[idx]

    barcode_umi <- data.table::data.table(
        barcode = barcode,
        taxid = as.integer(taxid),
        umi = umi
    )
    barcode_umi <- unique(barcode_umi)
    barcode_umi[, umi := NULL]
    barcode_umi[, umi := 1L]
    barcode_umi <- barcode_umi[
        , list(umi = sum(umi)),
        by = c("barcode", "taxid")
    ][order(-umi)]

    data.table::fwrite(barcode_umi,
        file = file_path(out_dir, sample, ext = "all.barcodes.txt"),
        sep = "\t", row.names = FALSE, col.names = TRUE
    )

    # create taxa counts data -----------------------------------
    ## use main id to match the umi data
    collected_data <- data.table::data.table(
        rank_name_pairs = mpa$V1,
        main_id = kr$V7
    )
    collected_data <- collected_data[main_id %in% barcode_umi$taxid] # nolint
    collected_data[, rank_name_pairs := strsplit( # nolint
        rank_name_pairs, "|",
        fixed = TRUE
    )]
    cli::cli_alert("Defining parent taxa")
    p2 <- progressr::progressor(
        steps = nrow(collected_data),
        auto_finish = TRUE
    )
    .tax_data. <- future.apply::future_lapply(
        collected_data$rank_name_pairs,
        function(x) {
            tax <- data.table::tstrsplit(x, split = "__", fixed = TRUE)
            data.table::setDT(tax)
            data.table::setnames(tax, c("rank", "name"))
            tax[, name := gsub("[^[:alnum:]]+", "_", name, perl = TRUE)] # nolint
            tax <- tax[c("k", "p", "c", "o", "f", "g", "s"), on = "rank"]
            n <- nrow(tax)
            for (i in seq_len(n)) {
                rank_i <- tax$rank[i]
                .value. <- c(
                    rep_len(NA_character_, i - 1L),
                    rep_len(tax$name[i], n - i + 1L)
                )
                tax[, c(rank_i) := .value.]
            }
            p2()
            tax
        }
    )
    collected_data[, tax_data := .tax_data.] # nolint

    collected_data[, rank_name_pairs := NULL] # nolint
    out <- collected_data[,
        data.table::rbindlist(tax_data, use.names = TRUE), # nolint
        by = "main_id"
    ][!is.na(name)] # nolint
    out <- merge(out, barcode_umi,
        by = c(main_id = "taxid"),
        allow.cartesian = TRUE
    )
    out[, main_id := NULL] # nolint
    out[, taxid := kr$V7[data.table::chmatch(name, kr$V8)]] # nolint
    data.table::setnames(
        out, c("k", "p", "c", "o", "f", "g", "s", "umi"),
        c("kingdom", "phylum", "class", "order", "family", "genus", "species", "counts")
    )
    out <- out[,
        list(counts = sum(counts, na.rm = TRUE)), # nolint
        by = c(
            "barcode", "taxid", "rank", "kingdom", "phylum", "class",
            "order", "family", "genus", "species"
        )
    ]
    data.table::setcolorder(
        out, c(
            "barcode", "taxid", "rank", "kingdom", "phylum", "class",
            "order", "family", "genus", "species"
        )
    )
    data.table::fwrite(out,
        file = file_path(out_dir, sample, ext = "counts.txt"),
        sep = "\t", row.names = FALSE, col.names = TRUE
    )
    invisible(out)
}
utils::globalVariables(c(
    "V8", "V1", "main_id", "rank_name_pairs", "name", "tax_data", "counts"
))
