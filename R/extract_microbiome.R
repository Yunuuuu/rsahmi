#' Extract Microbiome reads and Kraken output
#' @param fq1,fq2 The path to the classified microbiome fastq file (returned by
#' `run_kraken2`).
#' @param kraken_out The path to kraken output file, if `NULL`, will locate
#'   __sample.kraken.output.txt__ in `out_dir`.
#' @param kraken_report The path to kraken2uniq report, if `NULL`,
#'   will locate __sample.kraken.report.txt__ in `out_dir`.
#' @param mpa_report The path to mpa style report, if `NULL`,
#'   will locate __sample.kraken.report.mpa.txt__ in `out_dir`.
#' @param microbiome_pattern Perl regex patterns to match the bacteria
#'   microbiome.
#' @param ntaxid Number of taxids to extract at a time.
#' @inheritParams run_kraken2
#' @param ... Other arguments passed to [grepl].
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @export
extract_microbiome <- function(fq1, fq2 = NULL, kraken_out = NULL, kraken_report = NULL, mpa_report = NULL, out_dir = getwd(), sample = NULL, microbiome_pattern = "(?i)Bacteria|Fungi|Viruses", ntaxid = 8000L, ..., sys_args = list()) {
    sample <- sample %||% sub("_0*[12]?\\.(fastq|fq)(\\.gz)?$", "",
        basename(fq1),
        perl = TRUE
    )
    kraken_out <- define_path(kraken_out, sample = sample, out_dir = out_dir)
    kraken_report <- define_path(kraken_report,
        sample = sample,
        out_dir = out_dir
    )
    mpa_report <- define_path(mpa_report,
        sample = sample,
        out_dir = out_dir
    )
    if (!is.null(fq2)) {
        if (fq1 == fq2) {
            cli::cli_abort("{.arg fq1} and {.arg fq2} must be different.")
        }
    }
    taxid <- get_taxid(
        kraken_report = kraken_report,
        mpa_report = mpa_report,
        microbiome_pattern = microbiome_pattern,
        ...
    )
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    out_dir <- normalizePath(out_dir, mustWork = TRUE)
    microbiome_kraken_out(
        kraken_out = kraken_out,
        taxid = taxid, ntaxid = ntaxid,
        out_dir = out_dir,
        sample = sample,
        sys_args = sys_args
    )
    fq_files <- list(fq1 = fq1, fq2 = fq2)
    for (i in names(fq_files)) {
        fq <- fq_files[[i]]
        if (is.null(fq)) next
        cli::cli_alert("Extracting microbiome reads for {.arg {i}}")
        microbiome_reads(
            fq = fq, taxid = taxid,
            ntaxid = ntaxid,
            out_dir = out_dir,
            sys_args = sys_args
        )
    }
}
utils::globalVariables(c("r", "V1", "row_ids", "h"))

microbiome_kraken_out <- function(kraken_out, taxid, out_dir, sample = NULL, ntaxid = 8000L, sys_args = list()) {
    taxid_list <- split(taxid, ceiling(seq_along(taxid) / ntaxid))

    out_file <- file_path(out_dir, sample, ext = "microbiome.output.txt")
    if (file.exists(out_file)) file.remove(out_file)

    # extract microbiomme output -----------------------------------
    kraken_out <- normalizePath(kraken_out, mustWork = TRUE)
    out_file <- normalizePath(out_file, mustWork = FALSE)
    cli::cli_alert("Extracting microbiome kraken2 output")
    cli::cli_progress_bar(
        total = length(taxid_list),
        format = "{cli::pb_spin} Extracting kraken2 output | {cli::pb_current}/{cli::pb_total}",
        format_done = "Total time: {cli::pb_elapsed_clock}",
        clear = FALSE, auto_terminate = FALSE
    )
    for (i in seq_along(taxid_list)) {
        taxid <- paste0("(taxid ", taxid_list[[i]], ")", collapse = "\\|")
        run_command(
            args = c(
                sprintf("-w %s", shQuote(taxid)),
                kraken_out, ">>",
                out_file
            ),
            cmd = NULL,
            name = "grep",
            sys_args = sys_args,
            verbose = FALSE
        )
        cli::cli_progress_update()
    }
    cli::cli_process_done()
}

#' @param fq Path to the classified microbiome fastq file
#' @noRd
microbiome_reads <- function(fq, taxid, out_dir, ntaxid = 8000L, sys_args = list()) {
    taxid_list <- split(taxid, ceiling(seq_along(taxid) / ntaxid))
    line_number_file <- tempfile("line_number_", fileext = ".txt")
    if (file.exists(line_number_file)) {
        file.remove(line_number_file)
    }
    file.create(line_number_file)
    on.exit(file.remove(line_number_file))

    line_number_file <- normalizePath(line_number_file, mustWork = TRUE)
    fq <- normalizePath(fq, mustWork = TRUE)
    cli::cli_alert("Finding reads")
    cli::cli_progress_bar(
        total = length(taxid_list),
        format = "{cli::pb_spin} Finding reads | {cli::pb_current}/{cli::pb_total}",
        format_done = "Total time: {cli::pb_elapsed_clock}",
        clear = FALSE, auto_terminate = FALSE
    )
    for (i in seq_along(taxid_list)) {
        taxid <- paste0("taxid|", taxid_list[[i]], collapse = "\\|")
        run_command(
            args = c(
                "-wn", shQuote(taxid), fq, "|",
                "grep -Eo", shQuote("^[^:]+"), ">>",
                line_number_file
            ),
            cmd = NULL, name = "grep",
            sys_args = sys_args,
            verbose = FALSE
        )
        cli::cli_progress_update()
    }
    cli::cli_process_done()

    data <- data.table::fread(
        line_number_file,
        header = FALSE, sep = "\t",
        select = 1L
    )
    data[, r := V1 + 1L] # nolint
    data.table::setnames(data, c("h", "r"))
    data[, row_ids := seq_along(h)] # nolint
    data <- data.table::melt(
        data,
        id.vars = "row_ids",
        measure.vars = c("h", "r"),
        variable.name = "name",
        value.name = "value",
        variable.factor = FALSE
    )
    data.table::setorderv(data, c("row_ids", "name"), na.last = TRUE)
    data.table::fwrite(data[, "value"],
        file = line_number_file,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
    )

    cli::cli_alert("Extracting reads")
    out_file <- sub("\\.fq$", "", basename(fq), perl = TRUE)
    out_file <- file_path(out_dir, out_file, ext = "fa")
    run_command(
        c(
            shQuote("NR==FNR { lines[$1]; next } (FNR in lines)"),
            line_number_file, fq, ">", out_file
        ),
        cmd = NULL,
        name = "gawk"
    )
    run_command(
        c("-i", shQuote("s/@/>/g"), out_file),
        cmd = NULL, name = "sed"
    )
    cli::cli_alert_success("Extracting reads Done")
}

get_taxid <- function(kraken_report, mpa_report, microbiome_pattern, ...) {
    kr <- data.table::fread(
        kraken_report,
        sep = "\t",
        header = FALSE,
        select = 7L
    )[-c(1:2)]
    mpa <- data.table::fread(
        mpa_report,
        sep = "\t",
        header = FALSE,
        select = 1L
    )

    # get matched taxid
    taxid <- kr[[1L]][
        grepl(microbiome_pattern, mpa[[1L]], perl = TRUE, ...)
    ]
    taxid[!is.na(taxid)]
}
