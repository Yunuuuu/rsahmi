#' Sample level signal denoising
#'
#' In the low-microbiome biomass setting, real microbes also exhibit a
#' proportional number of total k-mers, number of unique k-mers, as well as
#' number of total assigned sequencing reads across samples; i.e. the following
#' three Spearman correlations are significant when tested using sample-level
#' data provided in Kraken reports: `cor(#k-mers, #unique k-mers)`,
#' `cor(#k-mers, #reads)` and `cor(#reads, #unique k-mers)`.
#'
#' @param kraken_reports The paths to kraken2uniq report for all samples, often
#'   in the form of __sample.kraken.report.txt__. It's easy to use [locate_path]
#'   to create the path.
#' @param samples An atomic character specifying sample names.
#' @inheritParams read_kr
#' @return A [data.table][data.table::data.table] of correlation coefficient and
#'   pvalue for `cor(#k-mers, #unique k-mers)` (r1 and p1), `cor(#k-mers,
#'   #reads)` (r2 and p2) and `cor(#reads, #unique k-mers)` (r3 and p3).
#' @export
sahmi_slsd <- function(kreports, method = "spearman", ...,
                       min_reads = 2L, min_uniq = 2L, min_number = 3L) {
    kreports <- pl$concat(kreports, how = "vertical")
    kreports[, .SD[.N > min_number], by = "taxid"][rank %in% c("G", "S")][,
        {
            min_vs_uniq <- stats::cor.test(
                min, uniq,
                method = method, ...
            )
            min_vs_reads <- stats::cor.test(
                min, reads,
                method = method, ...
            )
            reads_vs_uniq <- stats::cor.test(
                reads, uniq,
                method = method, ...
            )
            list(
                rank = unique(rank),
                r1 = min_vs_uniq$estimate,
                r2 = min_vs_reads$estimate,
                r3 = reads_vs_uniq$estimate,
                p1 = min_vs_uniq$p.value,
                p2 = min_vs_reads$p.value,
                p3 = reads_vs_uniq$p.value
            )
        },
        by = "name"
    ]
}
