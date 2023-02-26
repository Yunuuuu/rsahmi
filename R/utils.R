`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

run_command <- function(args = character(), cmd, name = NULL, sys_args = list(), verbose = TRUE) {
    if (!is.null(cmd)) {
        if (!file.exists(cmd) &&
            nzchar(Sys.which(cmd)) > 0L) {
            cmd <- Sys.which(cmd)
        }
        cmd <- normalizePath(cmd, mustWork = TRUE)
    } else if (!is.null(name)) {
        cmd <- Sys.which(name)
        if (nzchar(cmd) == 0L) {
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

file_path <- function(..., ext = NULL) {
    dots <- list(...)
    dots_len <- length(dots)
    if (dots_len > 0L && !is.null(ext)) {
        dots[[dots_len]] <- paste(dots[[dots_len]], ext, sep = ".")
    }
    do.call(file.path, dots)
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

column_to_rownames <- function(.data, var) {
    rownames(.data) <- as.character(.data[[var]])
    .data[[var]] <- NULL
    .data
}

str_match <- function(string, pattern, ignore.case = FALSE) {
    out <- regmatches(
        string,
        regexec(pattern, string,
            perl = TRUE, fixed = FALSE,
            ignore.case = ignore.case
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

str_extract <- function(string, pattern, ignore.case = FALSE) {
    matches <- regexpr(pattern, string,
        perl = TRUE, fixed = FALSE,
        ignore.case = ignore.case
    )
    start <- as.vector(matches)
    end <- start + attr(matches, "match.length") - 1L
    start[start == -1L] <- NA_integer_
    substr(string, start, end)
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
