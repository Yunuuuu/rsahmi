#' Refine Sequences from FASTQ Files with UMI/Barcode Actions
#'
#' This function refines one or two FASTQ files by applying trimming and
#' embedding operations defined via user-specified actions such as [`embed()`],
#' [`trim()`], or [`embed_trim()`]. It supports single-end and paired-end reads,
#' processes data in chunks, and uses multithreading for performance.
#'
#' @param fq1 Path to the first FASTQ file (required).
#' @param ofile1 Path to the output FASTQ file for `fq1` (required).
#' @param fq2 Optional path to the second (paired-end) FASTQ file. If you want
#'   to operate only on the second read, you may supply that file path via `fq1`
#'   instead.
#' @param ofile2 Optional path to the output FASTQ file for `fq2`.
#' @param umi_action1,umi_action2 Sequence action for extracting or trimming UMI
#' from `fq1`/`fq2`. `umi_action2` is only allowed if `fq2` is provided.
#' Specify sequence ranges using [`seq_range()`] or combine multiple ranges with
#' `c()`, this also applies to following all actions. By default, the UMI is
#' embedded in the header and trimmed from the sequence and quality
#' ([`embed_trim()`]). You can customize the behavior with [`embed()`], or
#' [`trim()`].
#' @param barcode_action1,barcode_action2 Sequence action for extracting or
#' trimming barcodes from `fq1`/`fq2`. `barcode_action2` is only allowed if
#' `fq2` is provided. By default, barcodes are embedded in the header and
#' trimmed from both the sequence and quality strings ([`embed_trim()`]).
#' Specify ranges as with UMI actions, and customize using [`embed()`], or
#' [`trim()`].
#' @param extra_actions1,extra_actions2 Additional sequence actions to apply to
#'   `fq1`/`fq2`. `extra_actions2` is only allowed if `fq2` is provided. These
#'   can be a single one or a list of them. By default, these actions
#'   perform trimming of sequences and qualities unless otherwise specified.
#' @param batch_size Integer. Number of records to accumulate before triggering
#' a write operation. Default is `r code_quote(BATCH_SIZE, quote = FALSE)`.
#' @param chunk_size Integer. Size in bytes of the intermediate chunk used to
#'   split and distribute data to worker threads during processing. Default is
#'   `10 * 1024 * 1024` (10MB).
#' @param buffer_size Integer specifying the buffer size in bytes used for
#' writing to disk. This controls the capacity of the buffered file writer.
#' Default is `1 * 1024 * 1024` (1MB).
#' @param nqueue Integer. Maximum number of buffers per thread, controlling the
#'   amount of in-flight data awaiting writing. Default: `3`.
#' @param threads Integer. Number of threads to use. Default will determined
#' atomatically by rayon.
#' @param odir A string of directory to save the output files. Please see
#' `Value` section for details.
#'
#' @return None. Outputs processed FASTQ files as specified by `ofile1` and
#' `ofile2`.
#' @details
#' Actions define what to do with sequence ranges specified using
#' [`seq_range()`].
#'
#' - UMI/barcode actions must be a single range (no lists) and default to
#'   [`embed_trim()`] if not specified.
#' - Extra actions can be multiple ranges (use a list) and default to [`trim()`]
#'   if not specified.
#' - Use [`embed()`], [`trim()`], or [`embed_trim()`] to specify the behavior.
#'
#' @export
seq_refine <- function(fq1, ofile1, fq2 = NULL, ofile2 = NULL,
                       umi_action1 = NULL, umi_action2 = NULL,
                       barcode_action1 = NULL, barcode_action2 = NULL,
                       extra_actions1 = NULL, extra_actions2 = NULL,
                       chunk_size = NULL, buffer_size = NULL,
                       batch_size = NULL, nqueue = NULL,
                       threads = NULL, odir = NULL,
                       mmap = TRUE) {
    assert_string(fq1, allow_empty = FALSE)
    assert_string(ofile1, allow_empty = FALSE, allow_null = TRUE)
    assert_string(fq2, allow_empty = FALSE, allow_null = TRUE)
    assert_string(ofile2, allow_empty = FALSE, allow_null = TRUE)
    umi_action1 <- check_ub_action(umi_action1, "UMI")
    barcode_action1 <- check_ub_action(barcode_action1, "BARCODE")
    extra_actions1 <- check_extra_actions(extra_actions1)

    umi_action2 <- check_ub_action(umi_action2, "UMI")
    barcode_action2 <- check_ub_action(barcode_action2, "BARCODE")
    extra_actions2 <- check_extra_actions(extra_actions2)
    if (is.null(fq2) &&
        (!is.null(ofile2) || !is.null(umi_action2) ||
            !is.null(barcode_action2) || !is.null(extra_actions2))) {
        cli::cli_abort(c(
            "The following arguments must not be used when {.arg fq2} is {.code NULL}:",
            "x" = "{.arg ofile2}",
            "x" = "{.arg umi_action2}",
            "x" = "{.arg barcode_action2}",
            "x" = "{.arg extra_actions2}",
            i = "These arguments are only applicable when paired-end reads are provided via {.arg fq2}."
        ))
    }

    if (is.null(ofile1) && is.null(ofile2)) {
        cli::cli_abort(c(
            "No output specified.",
            i = "Please provide at least one of {.arg ofile1} or {.arg ofile2} to write the results."
        ))
    }

    assert_number_whole(chunk_size, min = 1, allow_null = TRUE)
    assert_number_whole(buffer_size, min = 1, allow_null = TRUE)
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    assert_number_whole(threads,
        min = 1, max = as.double(parallel::detectCores()),
        allow_null = TRUE
    )
    nqueue <- check_queue(nqueue, 3L) *
        if (is.null(threads) || threads == 0L) {
            parallel::detectCores()
        } else {
            threads
        }
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    assert_bool(mmap)
    odir <- odir %||% getwd()
    dir_create(odir)
    chunk_size <- chunk_size %||% CHUNK_SIZE
    buffer_size <- buffer_size %||% BUFFER_SIZE
    batch_size <- batch_size %||% BATCH_SIZE
    threads <- threads %||% 0L # Let rayon to determine the threads
    actions1 <- c(list(umi_action1, barcode_action1), extra_actions1)
    actions1 <- actions1[
        !vapply(actions1, is.null, logical(1L), USE.NAMES = FALSE)
    ]
    if (length(actions1) == 0L) actions1 <- NULL
    actions2 <- c(list(umi_action2, barcode_action2), extra_actions2)
    actions2 <- actions2[
        !vapply(actions2, is.null, logical(1L), USE.NAMES = FALSE)
    ]
    if (length(actions2) == 0L) actions2 <- NULL
    if (is.null(actions1) && is.null(actions2)) {
        cli::cli_abort(c(
            "No sequence actions were specified.",
            i = "Please provide at least one action to proceed"
        ))
    }

    rust_call(
        "seq_refine",
        fq1 = fq1, ofile1 = ofile1,
        `_fq2` = fq2, `_ofile2` = ofile2,
        actions1 = actions1, `_actions2` = actions2,
        chunk_size = chunk_size,
        buffer_size = buffer_size,
        batch_size = batch_size,
        nqueue = nqueue,
        mmap = mmap,
        threads = threads
    )

    cli::cli_inform(c("v" = "Finished"))
}

check_ub_action <- function(action, tag, arg = caller_arg(action),
                            call = caller_env()) {
    if (is.null(action)) {
        action
    } else if (is_action(action)) {
        attr(action, "tag") <- attr(action, "tag", exact = TRUE) %||% tag
        action
    } else if (is_range(action)) {
        embed_trim(action, tag)
    } else {
        cli::cli_abort(c(
            "{.arg barcode_action2} must be created with {.fn seq_range} or a combination of them using {.fn c}.",
            i = "The default action is {.fn embed_trim}, but you can also use {.fn embed} or {.fn trim}."
        ), call = call)
    }
}

check_extra_actions <- function(actions, arg = caller_arg(actions),
                                call = caller_env()) {
    if (is.null(actions)) {
        return(actions)
    }

    if (is_range(actions)) {
        actions <- list(actions)
    } else if (is.list(actions)) {
        if (any(!vapply(actions, is_range, logical(1), USE.NAMES = FALSE))) {
            cli::cli_abort(c(
                "{.arg {arg}} must be created with {.fn seq_range} or a combination of them using {.fn c}.",
                i = "You can wrap multiple actions using {.fn list}.",
                i = "The default action is {.fn trim}, but you may also use {.fn embed} or {.fn embed_trim}."
            ), call = call)
        }
    } else {
        cli::cli_abort(c(
            "{.arg {arg}} must be created with {.fn seq_range} or a combination of them using {.fn c}.",
            i = "You can wrap multiple actions using {.fn list}.",
            i = "The default action is {.fn trim}, but you may also use {.fn embed} or {.fn embed_trim}."
        ), call = call)
    }
    lapply(actions, function(action) {
        if (!is_action(action)) action <- trim(action)
        action
    })
}
