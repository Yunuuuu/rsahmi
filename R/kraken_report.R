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
    kreport <- pl$scan_csv(kraken_report, separator = "\t", has_header = FALSE)$
        filter(pl$col("column_6")$neq("U"))$
        select(
        # rename necessary columns
        pl$col("column_7")$cast(pl$String)$alias("taxids"),
        pl$col("column_8")$str$strip_chars()$alias("taxon"),
        pl$col("column_6")$alias("ranks"),
        # kraken2 use the prefix blank space to specify the level
        pl$col("column_8")$str$extract(pl$lit("^( *)"), 1L)$
            str$len_chars()$div(2L)$alias("levels"),
        # rename necessary columns
        pl$col("column_1")$str$strip_chars()$cast(pl$Float64)$alias("percents"),
        pl$col("column_2")$alias("reads")$cast(pl$Int64)
    )$
        collect()

    kreport <- pl$DataFrame(
        parse_kreport_internal(
            kreport$height, kreport$to_list(), intermediate_ranks
        )
    )

    if (mpa) {
        kreport <- kreport$select(
            pl$col("taxon")$list$last()$alias("taxa"),
            pl$col("taxids")$list$last()$alias("taxid"),
            pl$col("ranks")$list$last()$alias("rank"),
            pl$col("percents")$list$last(),
            pl$col("reads")$list$last(),
            pl$col("ranks", "taxon")
        )$
            with_row_index("index")$
            explode(c("ranks", "taxon"))$
            # connect ranks with species names -------------
            group_by(
            c("index", "taxa", "taxid", "percents", "reads"),
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
            "taxa", "taxid",
            pl$col("phylogeny")$list$join("|"),
            "percents", "reads"
        )
    }
    kreport
}

# bench::mark(
#     dt = parse_kreport_data_table("kraken_report.txt"),
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
#     parse_kraken_report(kraken_report)$select(
#         pl$col("classification")$list$join("|")$alias("column_1"),
#         pl$col("reads")$alias("column_2")
#     )$to_data_frame(),
#     pl$read_csv(mpa_report, separator = "\t", has_header = FALSE)$
#         to_data_frame()
# )
# ✔ No differences
kraken_main_ranks <- c("R", "K", "D", "P", "C", "O", "F", "G", "S")
