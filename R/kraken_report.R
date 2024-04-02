#' Parse kraken report file
#'
#' @param kraken_report The path to kraken report file.
#' @param intermediate_ranks A bool indicates whether to include
#' non-traditional taxonomic ranks in output.
#' @param mpa A bool indicates whether to use mpa-style.
#' @return A polars [DataFrame][polars::DataFrame_class]
#' @importFrom polars pl
#' @export
parse_kraken_report <- function(kraken_report,
                                intermediate_ranks = TRUE,
                                mpa = FALSE) {
    # https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
    # 1. Percentage of fragments covered by the clade rooted at this taxon
    # 2. Number of fragments covered by the clade rooted at this taxon
    # 3. Number of fragments assigned directly to this taxon
    # * 4. Number of minimizers in read data associated with this taxon (new)
    # * 5. An estimate of the number of distinct minimizers in read data
    #    associated with this taxon (new)
    # 6. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
    #    (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that
    #    are not at any of these 10 ranks have a rank code that is formed by
    #    using the rank code of the closest ancestor rank with a number
    #    indicating the distance from that rank. E.g., "G2" is a rank code
    #    indicating a taxon is between genus and species and the grandparent
    #    taxon is at the genus rank.
    # 7. NCBI taxonomic ID number
    # 8. Indented scientific name
    kreport <- pl$scan_csv(kraken_report, separator = "\t", has_header = FALSE)$
        rename(
        percents = "column_1", total_reads = "column_2", reads = "column_3"
    )
    if (kreport$width == 6L) {
        cols <- pl$col("reads", "total_reads", "percents")
        kreport <- kreport$
            select(
            # rename necessary columns
            pl$col("column_5")$alias("taxids"),
            pl$col("column_6")$alias("taxon"),
            pl$col("column_4")$alias("ranks"),
            pl$col("reads"), pl$col("total_reads"), pl$col("percents")
        )
    } else if (kreport$width == 8L) {
        cols <- pl$col(
            "minimizer_len", "minimizer_n_unique",
            "reads", "total_reads", "percents"
        )
        kreport <- kreport$
            select(
            # rename necessary columns
            pl$col("column_7")$alias("taxids"),
            pl$col("column_8")$alias("taxon"),
            pl$col("column_6")$alias("ranks"),
            pl$col("column_4")$alias("minimizer_len"),
            pl$col("column_5")$alias("minimizer_n_unique"),
            pl$col("reads"), pl$col("total_reads"), pl$col("percents")
        )
    } else {
        cli::cli_abort("Invalid kraken2 report")
    }
    kreport <- kreport$
        filter(pl$col("ranks")$neq("U"))$
        with_columns(
        pl$col("percents")$str$strip_chars()$cast(pl$Float64),
        pl$col("taxids")$cast(pl$String),
        pl$col("taxon")$str$strip_chars(),
        # kraken2 use the prefix blank space to specify the level
        pl$col("taxon")$str$extract(pl$lit("^( *)"), 1L)$
            str$len_chars()$div(2L)$alias("levels")
    )$
        collect()
    kreport <- pl$DataFrame(
        parse_kreport_internal(
            kreport$height, kreport$to_list(), intermediate_ranks
        )
    )$
        with_columns(
        pl$col("reads")$cast(pl$List(pl$Int64)),
        pl$col("total_reads")$cast(pl$List(pl$Int64))
    )

    if (mpa) {
        kreport <- kreport$select(
            pl$col("taxon")$list$last()$alias("taxa"),
            pl$col("taxids")$list$last()$alias("taxid"),
            pl$col("ranks")$list$last()$alias("rank"),
            cols$list$last(), pl$col("ranks", "taxon")
        )$
            with_row_index("index")$
            explode(c("ranks", "taxon"))$
            # connect ranks with species names -------------
            group_by(
            c("index", "taxa", "taxid", "rank"),
            cols,
            maintain_order = TRUE
        )$
            agg(
            pl$concat_str(
                pl$when(pl$col("ranks")$is_in(kraken_main_ranks))$
                    then(pl$col("ranks")$str$to_lowercase())$
                    otherwise(pl$lit("x")),
                pl$col("taxon"),
                separator = "__"
            )$alias("phylogeny")
        )$
            # relocate columns -----------------------------
            select(
            "taxa", "taxid", "rank",
            pl$col("phylogeny")$list$join("|"), cols
        )
    }
    kreport
}

# bench::mark(
#     dt = parse_kraken_data_table("kraken_report.txt"),
#     polars = parse_kraken_report_polars("kraken_report.txt"),
#     check = FALSE
# )
# # A tibble: 2 × 13
#   expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc
#   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>
# 1 dt            334ms    365ms      2.74   183.5MB     8.22     2     6
# 2 polars        217ms    240ms      4.25     1.5MB     2.83     3     2
# parse_kreport_data_table <- function(kraken_report,
#                                      intermediate_ranks = TRUE) {
#     kr <- data.table::fread(
#         kraken_report,
#         sep = "\t", header = FALSE,
#         strip.white = FALSE
#     )
#     kr <- kr[V6 != "U", list(
#         ranks = V6,
#         taxids = as.character(V7),
#         phylogeny = paste(
#             ifelse(V6 %in% kraken_main_ranks, tolower(V6), "x"),
#             str_trim(V8),
#             sep = "__"
#         ),
#         levels = nchar(str_extract(V8, "^( *)")) / 2L,
#         percents = as.numeric(str_trim(V1)), reads = V2
#     )]
#     reports <- parse_kreport_internal(nrow(kr), kr, intermediate_ranks)
#     data.table::setDT(reports)[]
# }
# utils::globalVariables(c(paste0("V", c(1:2, 6:8))))

parse_kreport_internal <- function(n, kreport, intermediate_ranks) {
    # intermediate object to use
    prev_level <- NULL
    nms <- names(kreport)
    reports <- vector("list", n)
    report <- vector("list", length(kreport))
    names(report) <- nms
    for (i in seq_len(n)) {
        row <- lapply(kreport, .subset, i)
        # Move back ancestors if needed --------------
        actual_prev_level <- .subset2(row, "levels") - 1L
        if (!is.null(prev_level) && actual_prev_level != prev_level) {
            report <- lapply(
                report, .subset,
                .subset2(report, "levels") <= actual_prev_level
            )
        }
        # add current items into report --------------
        for (nm in nms) {
            report[[nm]] <- c(.subset2(report, nm), .subset2(row, nm))
        }
        if (!is.null(prev_level)) {
            # save current report --------------------
            if (intermediate_ranks ||
                any(.subset2(row, "ranks") == kraken_main_ranks)) {
                # don't save root species
                ranks <- .subset2(report, "ranks")
                keep <- ranks != "R"
                if (!intermediate_ranks) {
                    keep <- keep & ranks %in% kraken_main_ranks
                }
                reports[[i]] <- lapply(report, .subset, keep)
            }
        }
        prev_level <- .subset2(row, "levels")
    }
    reports <- reports[!vapply(reports, is.null, logical(1L))]
    reports <- lapply(nms, function(nm) lapply(reports, .subset2, nm))
    names(reports) <- nms
    reports
}

# parse_kraken_report("kraken_report.txt")
# waldo::compare(
#     rsahmi::parse_kraken_report(kraken_report, mpa = TRUE)$select(
#         pl$col("phylogeny")$alias("column_1"),
#         pl$col("total_reads")$alias("column_2")
#     )$to_data_frame(),
#     polars::pl$read_csv(mpa_report, separator = "\t", has_header = FALSE)$to_data_frame()
# )
# ✔ No differences
kraken_main_ranks <- c("R", "K", "D", "P", "C", "O", "F", "G", "S")
