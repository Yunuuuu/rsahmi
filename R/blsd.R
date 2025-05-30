#' Barcode level signal denoising
#'
#' True taxa are detected on multiple barcodes and with a proprotional number of
#' total and unique k-mer sequences across barcodes, measured as a significant
#' Spearman correlation between the number of total and unique k-mers across
#' barcodes. (`padj < 0.05`)
#'
#' @param kmer kmer data returned by [`prep_dataset()`].
#' @param method A character string indicating which correlation coefficient is
#'   to be used for the test. One of "pearson", "kendall", or "spearman", can be
#'   abbreviated.
#' @param ... Other arguments passed to [cor.test][stats::cor.test].
#' @param p.adjust Pvalue correction method, a character string. Can be
#'   abbreviated. Details see [p.adjust][stats::p.adjust].
#' @param min_kmer_len An integer, the minimal number of kmer to filter taxa.
#' SAHMI use `2`.
#' @param min_number An integer, the minimal number of cell per taxid. SAHMI use
#' `4`.
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @return A polars [DataFrame][polars::DataFrame_class]
#' @examples
#' \dontrun{
#' # 1. `sahmi_datasets` should be the output of all samples from 
#'       `prep_dataset()` 
#' # 2. `real_taxids_slsd` should be the output of `slsd()`
#' umi_list <- lapply(sahmi_datasets, function(dataset) {
#'     # Barcode level signal denoising (barcode k-mer correlation test)
#'     blsd <- blsd(dataset$kmer)
#'     real_taxids <- blsd$filter(pl$col("padj")$lt(0.05))$get_column("taxid")
#'     # only keep taxids pass Sample level signal denoising
#'     real_taxids <- real_taxids$filter(real_taxids$is_in(real_taxids_slsd))
#'     # remove contaminants
#'     real_taxids <- real_taxids$filter(
#'         real_taxids$is_in(attr(truly_microbe, "truly"))
#'     )
#'     # filter UMI data
#'     dataset$umi$filter(pl$col("taxid")$is_in(real_taxids))
#' })
#' }
#' @export
blsd <- function(kmer, method = "spearman", ..., p.adjust = "BH",
                 min_kmer_len = 3L, min_number = 3L) {
    use_polars()
    data_list <- kmer$
        filter(pl$col("kmer_len")$gt_eq(min_kmer_len))$
        filter(pl$len()$over("taxid", "taxa")$gt_eq(min_number))$
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
