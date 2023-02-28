#' Read kraken report
#' @param kraken_report The path to kraken2uniq report, often in the form of
#'   __sample.kraken.report.txt__.
#' @param sample A string, sample name.
#' @param min_reads Minimal number of reads.
#' @param min_uniq Minimal unique number.
#' @return A [data.table][data.table::data.table] object
#' @export 
read_kr <- function(kraken_report, sample = NULL, min_reads = 2L, min_uniq = 2L) {
    x <- data.table::fread(kraken_report, header = FALSE, sep = "\t")
    x[, V8 := str_trim(V8)] # nolint
    total_reads <- x$V2[1L] + x$V2[2L]
    n_microbiome_reads <- x[, sum(
        V2[V8 %in% c("Bacteria", "Fungi", "Viruses")], # nolint
        na.rm = TRUE
    )]
    out <- x[, list(
        sample = sample,
        rank = V6, taxid = V7, name = V8, # nolint
        reads = V2, min = V4, uniq = V5, # nolint
        rpm = V2 / total_reads * 1e6L,
        rpmm = V2 / n_microbiome_reads * 1e6L
    )]
    out[reads >= min_reads & uniq >= min_uniq] # nolint
}
utils::globalVariables(c(
    "V2", "V4", "V6", "V7", "reads", "uniq"
))
