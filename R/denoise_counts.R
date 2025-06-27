#' Denoise and Quantify Taxa from Kraken2 Results
#'
#' This function integrates three signal denoising strategies—cell-line-based
#' contaminant filtering, sample-level correlation filtering, and barcode-level
#' correlation testing— to identify reliable microbial taxa and quantify their
#' abundance across samples. It returns a table of taxa counts (per barcode),
#' after filtering out likely contaminants and false positives.
#'
#' @param kreports A character vector of paths to Kraken2 report files,
#' typically output by [`kuactor()`], one per sample.
#' @param kmers A character vector of paths to k-mer quantification files,
#' typically output by [`kuactor()`], one per sample.
#' @param umis A character vector of paths to UMI quantification files,
#' typically output by [`kuactor()`], one per sample.
#' @param samples A character of sample identifier for each element in `umis`.
#' If not provided, the names of the `umis` vector will be used.
#' @param cor_threshold Minimum correlation coefficient required in sample-level
#'   filtering. Default is `0`.
#' @param p_threshold Significance threshold for sample-level correlation
#'   p-values. Default is `0.05`.
#' @param padj_threshold Adjusted p-value threshold for barcode-level denoising.
#'   Default is `0.05`.
#' @inheritParams rpmm_quantile
#' @inheritParams slsd
#' @inheritParams blsd
#' @inheritParams taxa_counts
#' @export
denoise_counts <- function(kreports,
                           kmers,
                           umis, samples = NULL,
                           taxon = c("d__Bacteria", "d__Fungi", "d__Viruses"),
                           drop_unmatched_taxa = TRUE,
                           quantile = 0.95, alpha = 0.05,
                           min_reads = 3L,
                           min_minimizer_n_unique = 3L,
                           min_samples = 3L,
                           min_kmer_len = 3L,
                           min_cells = 3L,
                           cor_threshold = 0,
                           p_threshold = 0.05,
                           padj_threshold = 0.05,
                           ...,
                           study = "current study",
                           alternative = "greater",
                           p.adjust = "BH",
                           method = "spearman") {
    quantile_test_out <- rpmm_quantile(
        kreports = kreports,
        study = study,
        taxon = taxon,
        quantile = quantile,
        alpha = alpha,
        alternative = alternative,
        drop_unmatched_taxa = drop_unmatched_taxa
    )
    quantile_test_taxids <- attr(quantile_test_out, "truly", exact = TRUE)
    slsd_out <- slsd(
        kreports,
        method = method,
        min_reads = min_reads,
        min_minimizer_n_unique = min_minimizer_n_unique,
        min_samples = min_samples,
        ...
    )
    slsd_taxids <- slsd_out$filter(
        pl$col("r1")$gt(cor_threshold),
        pl$col("r2")$gt(cor_threshold),
        pl$col("r3")$gt(cor_threshold),
        pl$col("p1")$lt(p_threshold),
        pl$col("p2")$lt(p_threshold),
        pl$col("p3")$lt(p_threshold)
    )$get_column("taxid")
    umi_list <- .mapply(function(kmer, umi) {
        # Barcode level signal denoising (barcode k-mer correlation test)
        blsd_data <- blsd(kmer,
            method = method, ...,
            min_kmer_len = min_kmer_len,
            min_cells = min_cells,
            p.adjust = p.adjust
        )
        umi_data <- pl$read_parquet(umi)
        real_taxids <- blsd_data$
            filter(pl$col("padj")$lt(padj_threshold))$
            get_column("taxid")

        # only keep taxids pass Sample level signal denoising
        real_taxids <- real_taxids$filter(real_taxids$is_in(slsd_taxids))

        # remove contaminants
        real_taxids <- real_taxids$filter(
            real_taxids$is_in(quantile_test_taxids)
        )

        # filter UMI data
        umi_data$umi$filter(pl$col("taxid")$is_in(real_taxids))
    }, list(kmers, umis), NULL)
    counts <- taxa_counts(umi_list, samples %||% names(umis))
    structure(
        list(
            counts = counts, slsd = slsd_out,
            quantile_test = quantile_test_out
        ),
        class = "rsahmi"
    )
}
