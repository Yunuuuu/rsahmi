#' Quantify K-mers and UMI Data by Taxon
#'
#' This function counts total and unique k-mers per taxon across cell barcodes,
#' using both the cell barcode and unique molecular identifier (UMI) to resolve
#' read identity at the single-cell level. It aggregates k-mer counts for each
#' taxonomic rank of interest (by default, genus and species), including all
#' descendant taxa within those ranks.
#'
#' @param koutreads Path to the output file produced by [`koutreads()`].
#' @inheritParams koutreads
#' @param umi_tag (Optional) A string specifying the tag used to extract unique
#' molecular identifiers (UMIs) from each read. If `NULL`, all reads are counted
#' as total fragments.  Otherwise, only unique UMIs per (barcode, taxon) are
#' counted.
#' @param barcode_tag (Optional) A string specifying the tag used to extract the
#' cell barcode from each read. If `NULL`, all reads are assumed to originate
#' from a single cell.
#' @export
krcount <- function(koutreads, kreport,
                    umi_tag = NULL, barcode_tag = NULL,
                    taxonomy = c("D__Bacteria", "D__Fungi", "D__Viruses"),
                    batch_size = NULL,
                    nqueue = NULL) {
    rust_krcount(
        koutreads = koutreads, kreport = kreport,
        umi_tag = umi_tag, barcode_tag = barcode_tag,
        taxonomy = taxonomy, batch_size = batch_size,
        nqueue = nqueue
    )
}

rust_krcount <- function(koutreads, kreport,
                         umi_tag = NULL, barcode_tag = NULL,
                         taxonomy = c("D__Bacteria", "D__Fungi", "D__Viruses"),
                         batch_size = NULL,
                         nqueue = NULL, odir = NULL, pprof = NULL) {
    assert_string(koutreads, allow_empty = FALSE, allow_null = FALSE)
    assert_string(kreport, allow_empty = FALSE, allow_null = FALSE)
    assert_string(umi_tag, allow_empty = FALSE, allow_null = TRUE)
    assert_string(barcode_tag, allow_empty = FALSE, allow_null = TRUE)
    if (!is.null(taxonomy)) {
        taxonomy <- as.character(taxonomy)
        taxonomy <- taxonomy[!is.na(taxonomy)]
        if (length(taxonomy) == 0L) taxonomy <- NULL
    }
    assert_number_whole(batch_size, min = 1, allow_null = TRUE)
    nqueue <- check_queue(nqueue, 3L, 1)
    assert_string(pprof, allow_empty = FALSE, allow_null = TRUE)
    batch_size <- batch_size %||% KOUTPUT_BATCH

    if (is.null(pprof)) {
        rust_call(
            "krcount",
            koutreads = koutreads, kreport = kreport,
            umi_tag = umi_tag, barcode_tag = barcode_tag,
            taxonomy = taxonomy, batch_size = batch_size,
            nqueue = nqueue
        )
    } else {
        assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
        odir <- odir %||% getwd()
        rust_call(
            "pprof_krcount",
            koutreads = koutreads, kreport = kreport,
            umi_tag = umi_tag, barcode_tag = barcode_tag,
            taxonomy = taxonomy, batch_size = batch_size,
            nqueue = nqueue,
            pprof_file = file.path(odir, pprof)
        )
    }
}
