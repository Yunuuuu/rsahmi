#' Define sequence range behavior: embed, trim, or both
#'
#' These functions annotate a range or a list of ranges to indicate how
#' the specified subsequences (e.g. UMI or barcode) should be handled during
#' FASTQ processing.
#'
#' - `embed()` extracts the sequence and appends it to the read header without
#'   modifying the sequence.
#' - `trim()` removes the specified subsequence from the read but does not store
#'   it.
#' - `embed_trim()` both extracts the subsequence for the header and removes it
#'   from the read.
#'
#' Each function wraps the input range object in a new class to indicate its
#' behavior downstream.
#'
#' @param tag An character label used when embedding sequence content
#'   into the FASTQ header (used with `embed()` and `embed_trim()`).
#'
#'   For UMI and barcode actions, the tag will be assigned automatically as
#'   `"UMI"` and `"BARCODE"` respectively.
#'
#'   For other types of actions, you must explicitly specify a `tag` to ensure
#'   clarity in the embedded header.
#'
#' @param ranges A range or a list of ranges specifying the subsequence(s) to
#' process. Must be created using the [`seq_range()`] function.
#'
#' @return An annotated `rsahmi_seq_range` or `rsahmi_seq_ranges` object. object
#' with behavior-specific class:
#' \itemize{
#'   \item `rsahmi_embed` — for embedding only
#'   \item `rsahmi_trim` — for trimming only
#'   \item `rsahmi_embed_trim` — for embedding and trimming
#' }
#' @name subseq_actions
#' @export
embed <- function(tag, ranges) {
    assert_string(tag, allow_empty = FALSE, allow_null = FALSE)
    UseMethod("embed", ranges)
}

#' @export
embed.rsahmi_seq_range <- function(tag, ranges) {
    structure(
        ranges,
        tag = tag,
        class = c("rsahmi_embed", "rsahmi_seq_action", "rsahmi_seq_range")
    )
}

#' @export
embed.rsahmi_seq_ranges <- function(tag, ranges) {
    structure(
        ranges,
        tag = tag,
        class = c("rsahmi_embed", "rsahmi_seq_action", "rsahmi_seq_ranges")
    )
}

#' @rdname subseq_actions
#' @export
trim <- function(ranges) UseMethod("trim", ranges)

#' @export
trim.rsahmi_seq_range <- function(ranges) {
    structure(
        ranges,
        class = c("rsahmi_trim", "rsahmi_seq_action", "rsahmi_seq_range")
    )
}

#' @export
trim.rsahmi_seq_ranges <- function(ranges) {
    structure(
        ranges,
        class = c("rsahmi_trim", "rsahmi_seq_action", "rsahmi_seq_ranges")
    )
}

#' @rdname subseq_actions
#' @export
embed_trim <- function(tag, ranges) {
    assert_string(tag, allow_empty = FALSE, allow_null = FALSE)
    UseMethod("embed_trim", ranges)
}

#' @export
embed_trim.rsahmi_seq_range <- function(tag, ranges) {
    structure(
        ranges,
        tag = tag,
        class = c("rsahmi_embed_trim", "rsahmi_seq_action", "rsahmi_seq_range")
    )
}

#' @export
embed_trim.rsahmi_seq_ranges <- function(tag, ranges) {
    structure(
        ranges,
        tag = tag,
        class = c("rsahmi_embed_trim", "rsahmi_seq_action", "rsahmi_seq_ranges")
    )
}

#' @export
c.rsahmi_seq_action <- function(...) {
    cli::cli_abort(c(
        "Combining multiple {.cls rsahmi_seq_action} objects with {.fn c} is not supported.",
        i = "Use {.fn list} to collect multiple actions instead."
    ))
}

is_action <- function(action) inherits(action, "rsahmi_seq_action")
need_embed <- function(action) {
    inherits(action, c("rsahmi_embed_trim", "rsahmi_embed"))
}
