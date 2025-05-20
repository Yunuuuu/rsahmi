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
#' @param exclusive A boolean value, indicates whether taxa not found in
#' celllines data should be regarded as truly. Default: `FALSE`.
#' @return A polars [DataFrame][polars::DataFrame_class] with following
#' attributes:
#' 1. `pvalues`: Quantile test pvalue.
#' 2. `exclusive`: taxids in current study but not found in cellline data.
#' 3. `significant`: significant taxids with `pvalues < alpha`.
#' 4. `truly`: truly taxids based on `alpha` and `exclusive`. If `exclusive` is
#'    `TRUE`, this should be the union of `exclusive` and `significant`,
#'    otherwise, this should be the same with `significant`.
#' @examples
#' \dontrun{
#' # `paths` should be the output directory for each sample from
#' # `blit::kraken2()`
#' truly_microbe <- remove_contaminants(
#'     kraken_reports = file.path(paths, "kraken_report.txt"),
#'     quantile = 0.99, exclusive = FALSE
#' )
#' microbe_for_plot <- attr(truly_microbe, "truly")[
#'     order(attr(truly_microbe, "pvalue")[attr(truly_microbe, "truly")])
#' ]
#' microbe_for_plot <- microbe_for_plot[
#'     !microbe_for_plot %in% attr(truly_microbe, "exclusive")
#' ]
#' ggplot(
#'     truly_microbe$filter(pl$col("taxid")$is_in(microbe_for_plot))$
#'         to_data_frame(),
#'     aes(rpmm),
#' ) +
#'     geom_density(aes(fill = study), alpha = 0.5) +
#'     scale_x_log10() +
#'     facet_wrap(facets = vars(taxa), scales = "free") +
#'     theme(
#'         strip.clip = "off",
#'         axis.text = element_blank(),
#'         axis.ticks = element_blank(),
#'         legend.position = "inside",
#'         legend.position.inside = c(1, 0),
#'         legend.justification.inside = c(1, 0)
#'     )
#' }
#' @export
remove_contaminants <- function(kraken_reports, study = "current study",
                                taxon = c(
                                    "d__Bacteria", "d__Fungi", "d__Viruses"
                                ),
                                quantile = 0.95, alpha = 0.05,
                                alternative = "greater",
                                exclusive = FALSE) {
    use_polars()
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    cli::cli_alert_info("Parsing reads per million microbiome reads (rpmm)")
    kreports <- lapply(kraken_reports, parse_rpmm, taxon = taxon)
    kreports <- pl$concat(kreports, how = "vertical")

    # prepare celllines data ----------------------
    celllines <- pl$read_parquet(pkg_extdata("cell_lines.parquet"))$
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
    taxids <- vapply(rpmm_list, function(rpmm) {
        rpmm$slice(0L, 1L)$get_column("taxid")$to_r()
    }, character(1L))
    pvalues <- mapply(function(rpmm, taxid) {
        ref <- ref_quantile[taxid]
        if (is.na(ref)) return(NA_real_) # styler: off
        quantile_test(
            rpmm$get_column("rpmm")$to_r(),
            ref = ref,
            alternative = alternative
        )
    }, rpmm = rpmm_list, taxid = taxids, USE.NAMES = FALSE)
    exclusive_taxids <- setdiff(taxids, names(ref_quantile))
    truly <- significant <- taxids[!is.na(pvalues) & pvalues < alpha]
    if (isTRUE(exclusive)) {
        truly <- union(exclusive_taxids, truly)
    }

    # collect results and return ----------------
    structure(
        pl$concat(
            kreports$with_columns(study = pl$lit(study)),
            celllines$with_columns(study = pl$lit("cell lines")),
            how = "vertical"
        ),
        pvalues = structure(pvalues, names = taxids),
        exclusive = exclusive_taxids,
        significant = significant,
        truly = truly
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
