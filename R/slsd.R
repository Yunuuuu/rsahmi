#' Sample level signal denoising
#'
#' In the low-microbiome biomass setting, real microbes also exhibit a
#' proportional number of total k-mers, number of unique k-mers, as well as
#' number of total assigned sequencing reads across samples; i.e. the following
#' three Spearman correlations are significant when tested using sample-level
#' data provided in Kraken reports: `cor(minimizer_len, minimizer_n_unique)`,
#' `cor(minimizer_len, total_reads)` and `cor(total_reads, minimizer_n_unique)`.
#' (`r1>0 & r2>0 & r3>0 & p1<0.05 & p2<0.05 & p3<0.05`).
#'
#' @param kreports kreport files returned by [`kuactor()`] for all samples.
#' @param method A character string indicating which correlation coefficient is
#'   to be used for the test. One of "pearson", "kendall", or "spearman", can be
#'   abbreviated.
#' @param ... Other arguments passed to [cor.test][stats::cor.test].
#' @param min_reads An integer, the minimal number of the total reads to filter
#' taxa. SAHMI use `2`.
#' @param min_minimizer_n_unique An integer, the minimal number of the unique
#' number of minimizer to filter taxa. SAHMI use `2`.
#' @param min_samples An integer, the minimal number of samples per taxid. SAHMI
#' use `4`.
#' @return A polars [DataFrame][polars::DataFrame_class] of correlation
#'   coefficient and pvalue for `cor(minimizer_len, minimizer_n_unique)` (r1 and
#'   p1), `cor(minimizer_len, total_reads)` (r2 and p2) and `cor(total_reads,
#'   minimizer_n_unique)` (r3 and p3).
#' @export
slsd <- function(kreports, method = "spearman", ..., min_reads = 3L,
                 min_minimizer_n_unique = 3L, min_samples = 3L) {
    use_polars()
    kreports <- lapply(kreports, function(kreport) {
        pl$read_parquet(kreport)
    })
    data_list <- pl$concat(kreports, how = "vertical")$
        filter(
        pl$col("total_reads")$gt_eq(min_reads),
        pl$col("minimizer_n_unique")$gt_eq(min_minimizer_n_unique)
    )$
        filter(pl$len()$over("taxid", "taxa")$gt_eq(min_samples))$
        partition_by("taxid", "taxa")
    out_list <- lapply(data_list, function(data) {
        min_vs_uniq <- stats::cor.test(
            data$get_column("minimizer_len")$to_r(),
            data$get_column("minimizer_n_unique")$to_r(),
            method = method, ...
        )
        min_vs_reads <- stats::cor.test(
            data$get_column("minimizer_len")$to_r(),
            data$get_column("total_reads")$to_r(),
            method = method, ...
        )
        reads_vs_uniq <- stats::cor.test(
            data$get_column("total_reads")$to_r(),
            data$get_column("minimizer_n_unique")$to_r(),
            method = method, ...
        )
        data$select(pl$col("taxid", "taxa"))$slice(0L, 1L)$
            with_columns(
            r1 = min_vs_uniq$estimate,
            r2 = min_vs_reads$estimate,
            r3 = reads_vs_uniq$estimate,
            p1 = min_vs_uniq$p.value,
            p2 = min_vs_reads$p.value,
            p3 = reads_vs_uniq$p.value
        )
    })
    pl$concat(out_list, how = "vertical")
}
