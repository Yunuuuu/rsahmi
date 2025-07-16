#' Parse kraken report file
#'
#' @param kreport The path to kraken report file.
#' @param taxonomy A character vector. The set of taxonomic groups to include.
#' @return A data frame.
#' @seealso
#' <https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown>
#' @export
read_kreport <- function(kreport, taxonomy = NULL) {
    if (!is.null(taxonomy)) {
        taxonomy <- as.character(taxonomy)
        taxonomy <- taxonomy[!is.na(taxonomy)]
        if (length(taxonomy) == 0L) taxonomy <- NULL
    }
    out <- rust_call("read_kreport", kreport = kreport, taxonomy = taxonomy)
    class(out) <- "data.frame"
    attr(out, "row.names") <- .set_row_names(length(.subset2(out, 1L)))
    out
}
