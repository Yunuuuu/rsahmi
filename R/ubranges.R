#' Create a UMI/Barcode range
#'
#' Constructs a `rsahmi_ubrange` object representing a UMI/Barcode range with
#' optional start and end positions.
#'
#' @param start Integer. Start position (1-based). Optional if `end` is
#' provided.
#' @param end Integer. End position (1-based). Optional if `start` is provided.
#' @return A `rsahmi_ubrange` object.
#' @examples
#' # Create a single range
#' ubrange(start = 100, end = 200)
#'
#' # If `end` is NULL, the range goes to the end of the sequence
#' ubrange(start = 100)
#'
#' # If `start` is NULL, the range starts from the beginning of the sequence
#' ubrange(end = 200)
#'
#' # Combine multiple ranges
#' r1 <- ubrange(100, 200)
#' r2 <- ubrange(300, 400)
#' r3 <- ubrange(end = 500)
#'
#' # Combine them into a multi-range object using `c()`
#' c(r1, r2, r3)
#'
#' # Nested combinations are also supported
#' c(c(r1, r2), r3)
#' @export
ubrange <- function(start = NULL, end = NULL) {
    assert_number_whole(start, min = 1, allow_null = TRUE)
    assert_number_whole(end, min = 1, allow_null = TRUE)
    if (is.null(start) && is.null(end)) {
        cli::cli_abort("One of {.arg start} or {.arg end} must be provided")
    }
    new_ubrange(start, end)
}

new_ubrange <- function(start, end) {
    structure(list(start = start, end = end), class = "rsahmi_ubrange")
}

is_range <- function(x) inherits(x, c("rsahmi_ubrange", "rsahmi_ubranges"))

is_ubrange <- function(x) inherits(x, "rsahmi_ubrange")

is_ubranges <- function(x) inherits(x, "rsahmi_ubranges")

#' @export
c.rsahmi_ubrange <- function(...) {
    dots <- rlang::dots_list(..., .ignore_empty = "all")
    if (any(!vapply(dots, is_range, logical(1), USE.NAMES = FALSE))) {
        cli::cli_abort("{.cls ubrange} can only be combined with each other")
    }
    dots <- lapply(dots, function(e) {
        if (is_ubrange(e)) list(e) else e
    })
    structure(unlist(dots, FALSE, use.names = FALSE), class = "rsahmi_ubranges")
}

#' @export
c.rsahmi_ubranges <- c.rsahmi_ubrange

#' @export
print.rsahmi_ubrange <- function(x, ...) {
    cat("<ubrange>", "\n", sep = "")
    cat(
        .subset2(x, "start") %||% " ",
        .subset2(x, "end") %||% " ",
        sep = " -- "
    )
}

#' @export
print.rsahmi_ubranges <- function(x, ...) {
    cat("<ubranges>", "\n", sep = "")
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
