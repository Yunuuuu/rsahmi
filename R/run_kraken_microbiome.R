#' Run kraken2 and Extract Microbiome reads 
#' 
#' A helper function to run both `run_kraken2` and `extract_microbiome`
#' @param kraken2_cmd The path of kraken2 command.
#' @inheritParams extract_microbiome
#' @inheritParams run_kraken2
#' @export 
run_kraken_microbiome <- function(fq1, fq2 = NULL, sample = NULL, out_dir = getwd(), ncbi_blast_path = NULL, kraken2_db = NULL, cores = availableCores(), kraken2_cmd = NULL, kraken2_args = character(), python_cmd = NULL, microbiome_pattern = "(?i)Bacteria|Fungi|Viruses", ntaxid = 8000L, ..., sys_args = list()) {
    sample <- sample %||% sub("_0*[12]?\\.(fastq|fq)(\\.gz)?$", "",
        basename(fq1),
        perl = TRUE
    )
    kraken2_sys_args <- sys_args
    kraken2_sys_args$wait <- TRUE
    run_kraken2(
        fq1 = fq1, fq2 = fq2,
        sample = sample,
        out_dir = out_dir,
        ncbi_blast_path = ncbi_blast_path, # nolint
        kraken2_db = kraken2_db, # nolint
        mpa_report = TRUE,
        cores = cores,
        cmd = kraken2_cmd,
        kraken2_args = kraken2_args,
        python_cmd = python_cmd,
        sys_args = kraken2_sys_args
    )

    if (is.null(fq2)) {
        fq1 <- file_path(out_dir, sample, ext = "fq")
    } else {
        fq1 <- file_path(out_dir, sprintf("%s_1", sample), ext = "fq")
        fq2 <- file_path(out_dir, sprintf("%s_2", sample), ext = "fq")
    }
    extract_microbiome(
        fq1 = fq1, fq2 = fq2,
        kraken_out = file_path(out_dir, sample,
            ext = "kraken.output.txt"
        ),
        kraken_report = file_path(out_dir, sample,
            ext = "kraken.report.txt"
        ),
        mpa_report = file_path(out_dir, sample,
            ext = "kraken.report.mpa.txt"
        ),
        out_dir = out_dir,
        sample = sample, ...,
        microbiome_pattern = microbiome_pattern,
        ntaxid = ntaxid,
        sys_args = sys_args
    )
}
