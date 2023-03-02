#' Barcode level signal denoising
#'
#' True taxa are detected on multiple barcodes and with a proprotional number of
#' total and unique k-mer sequences across barcodes, measured as a significant
#' Spearman correlation between the number of total and unique k-mers across
#' barcodes.
#' @param kraken_report The path to kraken2uniq report, often in the form of
#'   __sample.kraken.report.txt__. It's easy to use [locate_path] to create the
#'   path. 
#' @param sckmer The path to sckmer results, often in the form of
#'   __sample.sckmer.txt__. It's easy to use [locate_path] to create the path. 
#' @param min_number An integer, the minimal number of items per taxid
#' @param method A character string indicating which correlation coefficient is
#'   to be used for the test. One of "pearson", "kendall", or "spearman", can be
#'   abbreviated.
#' @param ... Other arguments passed to [cor.test][stats::cor.test].
#' @param p.adjust Pvalue correction method, a character string. Can be
#'   abbreviated. Details see [p.adjust][stats::p.adjust].
#' @return A [data.table][data.table::data.table]
#' @export 
run_blsd <- function(kraken_report, sckmer, min_number = 3L, method = "spearman", ..., p.adjust = "BH") {
    kr_data <- data.table::fread(kraken_report, sep = "\t", header = FALSE)
    sckmer <- data.table::fread(sckmer, sep = "\t", header = TRUE)
    out <- sckmer[kmer > 1L, .SD[.N > min_number], by = "taxid"][, # nolint
        {
            cor_res <- stats::cor.test(
                kmer, uniq, # nolint
                method = method,
                ...
            )
            list(cor = cor_res$estimate, pvalue = cor_res$p.value)
        },
        by = "taxid"
    ]
    out[, padj := stats::p.adjust(pvalue, p.adjust)] # nolint
    out[, name := kr_data$V8[match(taxid, kr_data$V7)]]
    data.table::setcolorder(out, "name", after = "taxid")
    out[]
}
utils::globalVariables(c("kmer", "padj", "pvalue"))
