#' Create a UMI/Barcode range
#'
#' Constructs a `rsahmi_seq_range` object representing a UMI/Barcode range with
#' optional start and end positions.
#'
#' @param start Integer. Start position (1-based). Optional if `end` is
#' provided.
#' @param end Integer. End position (1-based). Optional if `start` is provided.
#' @return A `rsahmi_seq_range` or `rsahmi_seq_ranges` object.
#' @examples
#' # Create a single range
#' seq_range(start = 100, end = 200)
#'
#' # If `end` is NULL, the range goes to the end of the sequence
#' seq_range(start = 100)
#'
#' # If `start` is NULL, the range starts from the beginning of the sequence
#' seq_range(end = 200)
#'
#' # Combine multiple ranges
#' r1 <- seq_range(100, 200)
#' r2 <- seq_range(300, 400)
#' r3 <- seq_range(end = 500)
#'
#' # Combine them into a multi-range object using `c()`
#' c(r1, r2, r3)
#'
#' # Nested combinations are also supported
#' c(c(r1, r2), r3)
#' @export
seq_range <- function(start = NULL, end = NULL) {
    assert_number_whole(start, min = 1, allow_null = TRUE)
    assert_number_whole(end, min = 1, allow_null = TRUE)
    if (is.null(start) && is.null(end)) {
        cli::cli_abort("One of {.arg start} or {.arg end} must be provided")
    }
    if (!is.null(start)) start <- as.integer(start)
    if (!is.null(end)) end <- as.integer(end)
    new_seq_range(start, end)
}

new_seq_range <- function(start, end) {
    structure(
        list(start = start, end = end),
        class = "rsahmi_seq_range"
    )
}

is_range <- function(x) inherits(x, c("rsahmi_seq_range", "rsahmi_seq_ranges"))

is_seq_range <- function(x) inherits(x, "rsahmi_seq_range")

is_seq_ranges <- function(x) inherits(x, "rsahmi_seq_ranges")

#' @export
c.rsahmi_seq_range <- function(...) {
    dots <- rlang::dots_list(..., .ignore_empty = "all")
    if (any(!vapply(dots, is_range, logical(1), USE.NAMES = FALSE))) {
        cli::cli_abort("{.cls seq_range} can only be combined with each other")
    }
    dots <- lapply(dots, function(e) {
        if (is_seq_range(e)) list(e) else e
    })
    ranges <- unlist(dots, FALSE, use.names = FALSE)
    if (length(ranges) == 1L) {
        structure(.subset2(ranges, 1L), class = "rsahmi_seq_range")
    } else {
        structure(ranges, class = "rsahmi_seq_ranges")
    }
}

#' @export
c.rsahmi_seq_ranges <- c.rsahmi_seq_range

#' @export
print.rsahmi_seq_range <- function(x, ...) {
    cat("<seq_range>", "\n", sep = "")
    cat(
        .subset2(x, "start") %||% " ",
        .subset2(x, "end") %||% " ",
        sep = " -- "
    )
}

#' @export
print.rsahmi_seq_ranges <- function(x, ...) {
    cat(sprintf("<seq_ranges[%d]>", length(x)), "\n", sep = "")
    start <- vapply(x, function(range) {
        start <- .subset2(range, "start")
        if (is.null(start)) {
            ""
        } else {
            as.character(start)
        }
    }, character(1L), USE.NAMES = FALSE)
    end <- vapply(x, function(range) {
        end <- .subset2(range, "end")
        if (is.null(end)) {
            ""
        } else {
            as.character(end)
        }
    }, character(1L), USE.NAMES = FALSE)
    p <- paste0(
        format(c("start", start), justify = "right"),
        format(c("", rep_len(" -- ", length(start))), justify = "centre"),
        format(c("end", end), justify = "left")
    )
    cat(p, sep = "\n")
}
