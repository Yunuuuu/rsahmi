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

######################################################
column_to_rownames <- function(.data, var) {
    rownames(.data) <- as.character(.data[[var]])
    .data[[var]] <- NULL
    .data
}

handle_arg <- function(arg, name, format = "%s", sep = " ") {
    if (is.null(arg) || isFALSE(arg)) {
        return(NULL)
    } else if (isTRUE(arg)) {
        return(name)
    } else {
        return(sprintf(paste(name, format, sep = sep), arg))
    }
}

new_handlers <- function(message = "taxa processing") {
    progressr::handlers(progressr::handler_cli(
        format = sprintf(
            "{cli::pb_spin} %s | {cli::pb_current}/{cli::pb_total}", message
        ),
        format_done = "Total time: {cli::pb_elapsed_clock}",
        clear = FALSE
    ))
}

# file utils -----------------------------------------------
#' Locate the output file
#' @param x One of "kraken_report", "mpa_report", "kraken_out",
#'   "microbiome_out", "sckmer".
#' @param sample An atomic character specifying sample names, will be used to
#'   locate file since `SAHMI` use this to create output file name.
#' @param dir Path to retuls directory.
#' @return The path of x.
#' @export
locate_path <- function(x, sample, dir = getwd()) {
    locate_path_core(x, sample = sample, dir = dir)
}

define_path <- function(x, sample, dir = getwd()) {
    x %||% locate_path_core(deparse(substitute(x)),
        sample = sample, dir = dir
    )
}

locate_path_core <- function(x, sample, dir = getwd()) {
    x <- match.arg(x, c("kraken_report", "mpa_report", "kraken_out", "microbiome_out", "sckmer"))
    file_path(dir, sample, ext = switch_file(x))
}

switch_file <- function(x) {
    switch(x,
        kraken_report = "kraken.report.txt",
        mpa_report = "kraken.report.mpa.txt",
        kraken_out = "kraken.output.txt",
        microbiome_out = "microbiome.output.txt",
        sckmer = "sckmer.txt"
    )
}

file_path <- function(..., ext = NULL) {
    paths <- file.path(..., fsep = "/")
    if (!is.null(ext)) paths <- paste(paths, ext, sep = ".")
    paths
}
