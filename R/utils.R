`%||%` <- function(x, y) if (is.null(x)) y else x

run_command <- function(args = character(), cmd, name = NULL, sys_args = list(), verbose = TRUE) {
    if (!is.null(cmd)) {
        if (!file.exists(cmd) && nzchar(Sys.which(cmd))) {
            cmd <- Sys.which(cmd)
        }
        cmd <- normalizePath(cmd, mustWork = TRUE)
    } else if (!is.null(name)) {
        cmd <- Sys.which(name)
        if (!nzchar(cmd)) {
            cli::cli_abort("Cannot find {.field {name}} command")
        }
    }
    sys_args <- c(list(command = cmd, args = args), sys_args)
    if (verbose) {
        cli_args <- cli::cli_vec( # nolint
            args,
            list("vec-sep" = " ", "vec-last" = " ")
        )
        cli::cli_alert("Running command {.field {cmd} {cli_args}}")
    }
    do.call(system2, sys_args)
}

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

# string utils -----------------------------------------------
str_match <- function(string, pattern, ..., fixed = FALSE) {
    out <- regmatches(
        string,
        regexec(
            pattern = pattern, text = string,
            perl = !fixed, ..., fixed = fixed
        ),
        invert = FALSE
    )
    out <- lapply(out, function(x) {
        if (!length(x)) "" else x
    })
    out <- do.call("rbind", out)
    out[out == ""] <- NA_character_
    out
}

str_extract <- function(string, pattern, ..., fixed = FALSE) {
    matches <- regexpr(pattern, string, perl = !fixed, ..., fixed = fixed)
    start <- as.vector(matches)
    end <- start + attr(matches, "match.length") - 1L
    start[start == -1L] <- NA_integer_
    substr(string, start, end)
}

str_trim <- function(string, which = "both") {
    trimws(string, which = which, whitespace = "[\\h\\v]")
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
