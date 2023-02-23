extract_microbiome <- function(kraken_out, kraken_report, mpa_report, out_dir = getwd(), sample = NULL, microbiome_pattern = "(?i)Bacteria|Fungi|Viruses", ..., sys_args = list()) {
    taxid <- get_taxid(
        kraken_report = kraken_report,
        mpa_report = mpa_report,
        microbiome_pattern = microbiome_pattern,
        ...
    )
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    microbiome_out(
        kraken_out = kraken_out, taxid = taxid,
        out_dir = out_dir,
        sample = sample, sys_args = sys_args
    )
}

microbiome_out <- function(kraken_out, taxid, out_dir, sample = NULL, sys_args = list()) {
    taxid_list <- split(taxid, ceiling(seq_along(taxid) / 8000L))

    out_file <- file_path(out_dir, sample, ext = "microbiome.output.txt")
    if (file.exists(out_file)) file.remove(out_file)

    # extract microbiomme output -----------------------------------
    kraken_out <- normalizePath(kraken_out, mustWork = TRUE)
    out_file <- normalizePath(out_file, mustWork = FALSE)
    sys_args$wait <- TRUE
    cli::cli_alert("Extracting microbiome kraken2 output")
    cli::cli_progress_bar(
        total = length(taxid_list),
        format = "{cli::pb_spin} Extracting | {cli::pb_current}/{cli::pb_total}",
        format_done = "Total time: {cli::pb_elapsed_clock}"
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
            sys_args = sys_args
        )
        cli::cli_progress_update()
    }
    cli::cli_process_done()
}

microbiome_reads <- function() {
    
}

get_taxid <- function(kraken_report, mpa_report, microbiome_pattern, ...) {
    kr <- data.table::fread(
        kraken_report,
        sep = "\t",
        header = FALSE
    )
    kr <- kr[-c(1:2)]
    mpa <- data.table::fread(
        mpa_report,
        sep = "\t",
        header = FALSE
    )

    # get taxid
    taxid <- kr[[7L]][
        grepl(microbiome_pattern, mpa[[1L]], perl = TRUE, ...)
    ]
    taxid[!is.na(taxid)]
}
