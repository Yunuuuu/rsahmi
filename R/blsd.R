#' Barcode level signal denoising
#'
#' True taxa are detected on multiple barcodes and with a proprotional number of
#' total and unique k-mer sequences across barcodes, measured as a significant
#' Spearman correlation between the number of total and unique k-mers across
#' barcodes. (`padj < 0.05`)
#'
#' @param kmer kmer file returned by [`sckmer()`].
#' @param method A character string indicating which correlation coefficient is
#'   to be used for the test. One of "pearson", "kendall", or "spearman", can be
#'   abbreviated.
#' @param ... Other arguments passed to [cor.test][stats::cor.test].
#' @param min_kmer_len An integer, the minimal number of kmer to filter taxa.
#' SAHMI use `2`.
#' @param min_cells An integer, the minimal number of cell per taxid. SAHMI use
#' `4`.
#' @param p.adjust Pvalue correction method, a character string. Can be
#'   abbreviated. Details see [p.adjust][stats::p.adjust].
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @return A polars [DataFrame][polars::DataFrame_class]
#' @export
blsd <- function(kmer, method = "spearman", ...,
                 min_kmer_len = 3L, min_cells = 3L, p.adjust = "BH") {
    use_polars()
    kmer <- pl$read_parquet(kmer)
    data_list <- kmer$
        filter(pl$col("kmer_len")$gt_eq(min_kmer_len))$
        filter(pl$len()$over("taxid", "taxa")$gt_eq(min_cells))$
        partition_by("taxid", "taxa")
    out_list <- lapply(data_list, function(data) {
        cor_res <- stats::cor.test(
            x = data$get_column("kmer_len")$to_r(),
            y = data$get_column("kmer_n_unique")$to_r(), # nolint
            method = method,
            ...
        )
        data$select(pl$col("taxid", "taxa"))$slice(0L, 1L)$
            with_columns(cor = cor_res$estimate, pvalue = cor_res$p.value)
    })
    out <- pl$concat(out_list, how = "vertical")
    padj <- stats::p.adjust(out$get_column("pvalue")$to_r(), p.adjust)
    out$with_columns(padj = padj)
}
