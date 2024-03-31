`%||%` <- function(x, y) if (is.null(x)) y else x

is_scalar <- function(x) length(x) == 1L

# mimic polars list methods --------------------------
list_gather <- function(x, index, USE.NAMES = FALSE) {
    mapply(.subset, x = x, i = index, USE.NAMES = USE.NAMES, SIMPLIFY = FALSE)
}

list_get <- function(x, index, USE.NAMES = FALSE) {
    mapply(.subset, x = x, i = index, USE.NAMES = USE.NAMES)
}

list_last <- function(x, USE.NAMES = FALSE) {
    mapply(.subset, x = x, i = lengths(x), USE.NAMES = USE.NAMES)
}

list_first <- function(x, USE.NAMES = FALSE) {
    mapply(.subset, x = x, i = rep_len(1L, length(x)), USE.NAMES = USE.NAMES)
}

list_contains <- function(x, items) {
    vapply(x, function(xx) any(xx %in% items), logical(1L))
}
