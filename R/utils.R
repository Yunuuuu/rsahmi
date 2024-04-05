`%||%` <- function(x, y) if (is.null(x)) y else x

pkg_nm <- function() {
    utils::packageName(topenv(environment()))
}

internal_file <- function(...) {
    system.file(..., package = pkg_nm(), mustWork = TRUE)
}

is_scalar <- function(x) length(x) == 1L

dir_create <- function(path, ...) {
    if (!dir.exists(path) &&
        !dir.create(path = path, showWarnings = FALSE, ...)) {
        cli::cli_abort("Cannot create directory {.path {path}}")
    }
}

# purrr-like function -------------------------------
map2 <- function(.x, .y, .f, ..., USE.NAMES = FALSE) {
    mapply(.f, .x, .y, MoreArgs = list(...), # styler: off
        SIMPLIFY = FALSE, USE.NAMES = USE.NAMES
    )
}

imap <- function(.x, .f, ..., USE.NAMES = FALSE) {
    map2(.x, names(.x) %||% seq_along(.x), .f, ..., USE.NAMES = USE.NAMES)
}

transpose <- function(.l) {
    if (!length(.l)) return(.l) # styler: off

    inner_names <- names(.l[[1]])
    if (is.null(inner_names)) {
        fields <- seq_along(.l[[1]])
    } else {
        fields <- inner_names
        names(fields) <- fields
        .l <- lapply(.l, function(x) {
            if (is.null(names(x))) names(x) <- inner_names # styler: off
            x
        })
    }

    # This way missing fields are subsetted as `NULL` instead of causing
    # an error
    .l <- lapply(.l, as.list)

    lapply(fields, function(i) lapply(.l, .subset2, i))
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
