#' Identifying contaminants and false positives (cell line quantile test)
#'
#' @inheritParams extractor
#' @param quantile Probabilities with values in `[0, 1]` specifying the quantile
#' to calculate.
#' @return A polars [DataFrame][polars::DataFrame_class]
#' @export
contaminants <- function(kraken_report,
                         taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                         quantile = 0.99) {
    kreport <- parse_kraken_report(kraken_report)
    ref_reads <- kreport$
        filter(
        pl$concat_str(
            pl$col("ranks")$list$last()$str$to_lowercase(),
            pl$col("taxon")$list$last(),
            separator = "__"
        )$is_in(pl$lit(taxon))
    )$
        select(pl$col("total_reads")$sum())$to_series()

    rpmm <- kreport$
        with_row_index("index")$
        with_columns(
        pl$col("taxids")$list$last()$alias("taxid"),
        pl$col("total_reads")$list$last(),
        pl$col("taxon")$list$last()$alias("taxa")
    )$
        explode(pl$col("ranks", "taxon"))$
        filter(
        pl$concat_str(
            pl$col("ranks")$str$to_lowercase(),
            pl$col("taxon"),
            separator = "__"
        )$is_in(pl$lit(taxon))
    )$
        select("index", "taxid", "taxa", "total_reads")$
        unique()$
        select(
        pl$col("taxid", "taxa"),
        pl$col("total_reads")$div(ref_reads)$mul(10^6L)$alias("rpmm")
    )

    cell_lines <- pl$read_parquet("inst/extdata/cell_lines.parquet")$
        group_by("taxid")$
        agg(pl$col("rpmm")$log10()$quantile(quantile))$
        select(
        pl$col("taxid"),
        pl$col("rpmm")$alias("cell_lines_rpmm"),
        pl$lit(10)$pow(pl$col("rpmm"))$alias("cell_lines_rpmm_quantile")
    )
    rpmm$join(cell_lines, "taxid", how = "left")
}
