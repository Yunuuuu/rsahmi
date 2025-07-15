#' Benchmark File Reading with Optional Memory Mapping
#'
#' This function benchmarks the reading performance of a file by processing it
#' in chunks. It allows for optional memory-mapped file access, which can
#' improve performance in multi-threaded environments by reducing memory
#' copying.
#'
#' @inheritParams kractor
#' @export
bench_read <- function(file, chunk_size = NULL, mmap = TRUE) {
    assert_string(file, allow_empty = FALSE, allow_null = FALSE)
    assert_number_whole(chunk_size, min = 1, allow_null = TRUE)
    assert_bool(mmap)
    chunk_size <- chunk_size %||% CHUNK_SIZE
    rust_call("bench_read", file = file, chunk_size = chunk_size, mmap = mmap)
    cli::cli_inform(c("v" = "Finished"))
}

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

RUST_CALL <- .Call

#' @keywords internal
rust_method <- function(class, method, ...) {
    rust_call(sprintf("%s__%s", class, method), ...)
}

#' @keywords internal
rust_call <- function(.NAME, ..., call = caller_env()) {
    # call the function
    out <- RUST_CALL(sprintf("wrap__%s", .NAME), ...)

    # propagate error from rust --------------------
    if (!inherits(out, "extendr_result")) return(out) # styler: off
    if (!is.null(err <- .subset2(out, "err"))) {
        rlang::abort(err, call = call)
    }
    .subset2(out, "ok")
}
