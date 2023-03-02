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
#' @inheritParams run_blsd
#' @return A [data.table][data.table::data.table] of correlation coefficient and
#'   pvalue for `cor(#k-mers, #unique k-mers)` (r1 and p1), `cor(#k-mers,
#'   #reads)` (r2 and p2) and `cor(#reads, #unique k-mers)` (r3 and p3).
#' @export 
run_slsd <- function(kraken_reports, samples = NULL, method = "spearman", ..., min_number = 3L, min_reads = 2L, min_uniq = 2L) {
    if (is.null(samples)) {
        samples <- rep_len(list(NULL), length(kraken_reports))
    } else if (length(samples) != length(kraken_reports)) {
        cli::cli_abort("{.arg kraken_reports} and {.arg samples} must have the same length.")
    }
    kr_data <- .mapply(
        function(kraken_report, sample, ...) {
            read_kr(kraken_report, sample, ...)
        },
        list(kraken_report = kraken_reports, sample = samples),
        list(min_reads = min_reads, min_uniq = min_uniq)
    )
    kr_data <- data.table::rbindlist(kr_data, use.names = TRUE)

    kr_data[, .SD[.N > min_number], by = "taxid"][rank %in% c("G", "S")][,
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
