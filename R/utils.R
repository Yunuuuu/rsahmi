# mimic polars str methods ---------------------------
# https://rpolars.github.io/man/ExprStr_contains_any.html
str_contains_any <- function(string, patterns, ...) {
    str_detect(string = string, pattern = paste0(patterns, collapse = "|"))
}

is_scalar <- function(x) length(x) == 1L

dir_create <- function(path, ...) {
    if (!dir.exists(path) &&
        !dir.create(path = path, showWarnings = FALSE, ...)) {
        cli::cli_abort("Cannot create directory {.path {path}}")
    }
}

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
