#' Quantitation of microbes
#' @description After identifying true taxa, reads assigned to those taxa are
#' extracted and then undergo a series of filters. The cell barcode and UMI are
#' used to demultiplex the reads and create a barcode x taxa counts matrix. The
#' full taxonomic classification of all resulting barcodes and the number of
#' counts assigned to each clade are tabulated.
#'
#' @param umi_list A list of polars [DataFrame][polars::DataFrame_class] for UMI
#' data returned by [prep_dataset].
#' @param samples A character of sample identifier for each element in
#' `umi_list`.
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @export
taxa_counts <- function(umi_list, samples = NULL) {
    columns <- c("sample", "barcode", "taxid", "taxa", "rank")
    if (is.list(umi_list)) {
        if (!is.null(samples) && length(samples) != length(umi_list)) {
            cli::cli_abort(
                "{.arg samples} must have the same length of {.arg umi_list}"
            )
        }
        samples <- samples %||% names(umi_list) %||% seq_along(umi_list)
        umi_list <- .mapply(function(umi, sample) {
            umi$with_columns(sample = pl$lit(sample))
        }, list(umi = umi_list, sample = samples), NULL)
        umi <- pl$concat(umi_list, how = "vertical")
    } else if (polars::is_polars_df(umi_list)) {
        if (!is.null(samples)) {
            if (length(samples) != 1L) {
                cli::cli_abort(
                    "{.arg samples} must be a scalar string for a single umi"
                )
            }
            umi <- umi_list$with_columns(sample = pl$lit(samples))
        } else {
            columns <- setdiff(columns, "sample")
            umi <- umi_list
        }
    } else {
        cli::cli_abort(paste(
            "{.arg umi_list} must be a polars {.cls DataFrame}",
            "or a list of polars {.cls DataFrame}"
        ))
    }

    # create barcode umi data -----------------------------------
    counts <- umi$lazy()$
        with_columns(pl$col(pl$List(pl$String))$list$join("|"))$
        group_by(c(columns, "taxon", "ranks"))$
        agg(pl$col("umi")$n_unique())$
        with_columns(pl$col("ranks", "taxon")$str$split("|"))$
        #     with_columns(
        #     pl$col("ranks", "taxon")$list$eval(
        #         pl$int_range(0L, pl$element()$len())$add(1L)
        #     )$name$suffix("_len")
        # )$
        #     explode("^.+_len$")$
        #     with_columns(
        #     pl$col("ranks")$list$slice(0L, length = pl$col("ranks_len")),
        #     pl$col("taxon")$list$slice(0L, length = pl$col("taxon_len"))
        # )$
        #     drop("^.+_len$")$
        with_row_index("index")$
        explode(pl$col("ranks", "taxon"))$
        with_columns(
        # A rank code, indicating (U)nclassified, (R)oot, (D)omain,
        # (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or
        # (S)pecies. Taxa that are not at any of these 10 ranks have a rank
        # code that is formed by using the rank code of the closest ancestor
        # rank with a number indicating the distance from that rank. E.g.,
        # "G2" is a rank code indicating a taxon is between genus and
        # species and the grandparent taxon is at the genus rank.
        pl$col("ranks")$
            str$replace("^R", "root")$
            str$replace("^D", "domain")$
            str$replace("^K", "kingdom")$
            str$replace("^P", "phylum")$
            str$replace("^C", "class")$
            str$replace("^O", "order")$
            str$replace("^F", "family")$
            str$replace("^G", "genus")$
            str$replace("^S", "species")
    )$
        collect()$
        pivot(
        values = "taxon", index = c("index", columns, "umi"),
        columns = "ranks"
    )$
        drop("index")

    # relocate columns ----------------------------
    columns <- list(pl$col(columns))
    columns <- counts$columns
    for (taxa in c("root", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")) {
        pattern <- sprintf("^%s\\d*$", taxa)
        if (any(grepl(pattern, columns, perl = TRUE))) {
            columns <- c(columns, list(pl$col(pattern)))
        }
    }
    counts$select(c(columns, list(pl$col("umi")$alias("counts"))))
}
