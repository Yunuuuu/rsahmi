`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

run_command <- function(args = character(), cmd, name = NULL, sys_args = list()) {
    if (!is.null(cmd)) {
        if (!file.exists(cmd) &&
            nzchar(Sys.which(cmd)) > 0L) {
            cmd <- Sys.which(cmd)
        }
        cmd <- normalizePath(cmd, mustWork = TRUE)
    } else if (!is.null(name)) {
        cmd <- Sys.which(name)
        if (nzchar(cmd) == 0L) {
            cli::cli_abort("Cannot find {.fun {name}}")
        }
    }
    sys_args <- c(list(command = cmd, args = args), sys_args)
    cli_args <- cli::cli_vec( # nolint
        args,
        list("vec-trunc" = 3L, "vec-sep" = " ", "vec-last" = " ")
    )
    cli::cli_alert("Running command {.field {cmd} {cli_args}}")
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
