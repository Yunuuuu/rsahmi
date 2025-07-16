#' Parse kraken report file
#'
#' @param kreport The path to kraken report file.
#' @return A data frame.
#' @seealso
#' <https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown>
#' @export
read_kreport <- function(kreport) {
    out <- rust_call("read_kreport", kreport = kreport)
    class(out) <- "data.frame"
    attr(out, "row.names") <- .set_row_names(length(.subset2(out, 1L)))
    out
}
