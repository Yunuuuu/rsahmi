#' Identifying contaminants and false positives taxa (cell line quantile test)
#'
#' @param kraken_reports A character of path to all kraken report files.
#' @param study A string of the study name, used to differentiate with cell line
#' data.
#' @inheritParams extractor
#' @param quantile Probabilities with values in `[0, 1]` specifying the quantile
#' to calculate.
#' @param alpha Level of significance.
#' @param alternative A string specifying the alternative hypothesis, must be
#' one of "two.sided", "greater" (default) or "less". You can specify just the
#' initial letter.
#' @return A polars [DataFrame][polars::DataFrame_class] with following
#' attributes:
#' 1. `pvalues`: Quantile test pvalue.
#' 2. `contaminants`: significant taxids based on `alpha`.
#' @export
contaminants <- function(kraken_reports, study = "current study",
                         taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                         quantile = 0.95, alpha = 0.05,
                         alternative = "greater") {
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    cli::cli_alert_info("Parsing reads per million microbiome reads (rpmm)")
    kreports <- lapply(kraken_reports, parse_rpmm, taxon = taxon)
    kreports <- pl$concat(kreports, how = "vertical")

    # prepare celllines data ----------------------
    celllines <- pl$read_parquet(
        internal_file("extdata", "cell_lines.parquet")
    )$
        select(pl$col("taxid", "rpmm"))
    celllines <- kreports$select(pl$col("taxid", "taxa"))$
        join(celllines, on = "taxid", how = "inner")

    # Do quantile test ----------------------------
    cli::cli_alert_info("Doing quantile test")
    ref_quantile <- celllines$group_by("taxid")$
        agg(pl$col("rpmm")$quantile(quantile))$
        to_data_frame()
    ref_quantile <- structure(ref_quantile$rpmm, names = ref_quantile$taxid)
    rpmm_list <- kreports$partition_by("taxid", "taxa")
    names(rpmm_list) <- vapply(
        rpmm_list,
        function(x) x$slice(0L, 1L)$get_column("taxid")$to_r(),
        character(1L)
    )
    pvalues <- imap(rpmm_list, function(rpmm, taxid) {
        quantile_test(rpmm$get_column("rpmm")$to_r(),
            ref = ref_quantile[taxid],
            alternative = alternative
        )
    }, USE.NAMES = TRUE)

    # collect results and return ----------------
    structure(
        pl$concat(
            kreports$with_columns(study = pl$lit(study)),
            celllines$with_columns(study = pl$lit("cell lines")),
            how = "vertical"
        ),
        pvalues = pvalues,
        contaminants = names(pvalues)[!is.na(pvalues) & pvalues < alpha]
    )
}

parse_rpmm <- function(kraken_report, taxon) {
    kreport <- parse_kraken_report(kraken_report)
    ref_reads <- kreport$
        filter(
        pl$concat_str(
            pl$col("ranks")$list$last()$str$to_lowercase(),
            pl$col("taxon")$list$last(),
            separator = "__"
        )$is_in(pl$lit(taxon))
    )$
        select(pl$col("total_reads")$list$last()$sum())$
        to_series()$cast(pl$Float64)
    kreport$
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
        pl$col("total_reads")$alias("rpmm")$div(ref_reads)$mul(10^6L)
    )
}

# https://people.stat.sc.edu/hitchcock/Rexamples518section3_2.txt
# https://people.stat.sc.edu/hitchcock/notes518fall13sec32filledin.pdf
quantile_test <- function(x, ref = 0, p = .5, alternative) {
    n <- length(x)
    T1 <- sum(x <= ref)
    T2 <- sum(x < ref)
    switch(alternative,
        less = stats::pbinom(T2 - 1L, n, p, lower.tail = FALSE),
        greater = stats::pbinom(T1, n, p),
        two.sided = 2 * min(
            stats::pbinom(T2 - 1L, n, p, lower.tail = FALSE),
            stats::pbinom(T1, n, p)
        )
    )
}
