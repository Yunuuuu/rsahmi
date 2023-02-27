#' Run kraken2 and Extract Microbiome reads
#'
#' A helper function to run `run_kraken2`, `extract_microbiome` and `run_sckmer`
#' sequentially.
#' @param kraken2_cmd The path of kraken2 command.
#' @inheritParams extract_microbiome
#' @inheritParams run_kraken2
#' @inheritParams run_sckmer
#' @export
run_krmikm <- function(fq1, fq2 = NULL, sample = NULL, out_dir = getwd(),
                       ncbi_blast_path = NULL, kraken2_db = NULL,
                       cores = availableCores(), kraken2_cmd = NULL,
                       kraken2_args = character(), python_cmd = NULL,
                       microbiome_pattern = "(?i)Bacteria|Fungi|Viruses", ntaxid = 8000L, ranks = c("G", "S"),
                       cb_len = 16L, umi_len = 10L,
                       host = 9606L, nsample = 1000L,
                       kmer_len = 35L, min_frac = 0.5,
                       ..., sys_args = list()) {
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
        sample = sample,
        microbiome_pattern = microbiome_pattern,
        ntaxid = ntaxid, ...,
        sys_args = sys_args
    )
    fa1 <- sub("fq$", "fa", fq1, perl = TRUE)
    if (!is.null(fq2)) {
        fa2 <- sub("fq$", "fa", fq2, perl = TRUE)
    } else {
        fa2 <- fq2
    }
    run_sckmer(
        fa1 = fa1, fa2 = fa2,
        kraken_report = file_path(out_dir, sample,
            ext = "kraken.report.txt"
        ),
        mpa_report = file_path(out_dir, sample,
            ext = "kraken.report.mpa.txt"
        ),
        microbiome_output = file_path(out_dir, sample,
            ext = "microbiome.output.txt"
        ),
        sample = sample, out_dir = out_dir,
        ranks = ranks, cb_len = cb_len, umi_len = umi_len,
        host = host, nsample = nsample,
        kmer_len = kmer_len, min_frac = min_frac,
        cores = cores
    )
    invisible(list(
        fa1 = fa1, fa2 = fa2,
        kraken_report = file_path(out_dir, sample,
            ext = "kraken.report.txt"
        ),
        mpa_report = file_path(out_dir, sample,
            ext = "kraken.report.mpa.txt"
        ),
        kraken_out = file_path(out_dir, sample,
            ext = "kraken.output.txt"
        ),
        microbiome_output = file_path(out_dir, sample,
            ext = "microbiome.output.txt"
        ),
        kmer = file_path(out_dir, sample, ext = "sckmer.txt")
    ))
}
