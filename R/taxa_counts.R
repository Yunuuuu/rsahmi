#' Quantitation of microbes
#' @description After identifying true taxa, reads assigned to those taxa are
#' extracted and then undergo a series of filters. The cell barcode and UMI are
#' used to demultiplex the reads and create a barcode x taxa counts matrix. The
#' full taxonomic classification of all resulting barcodes and the number of
#' counts assigned to each clade are tabulated.
#'
#' @param umi UMI data returned by [sahmi_kmer].
#' @param taxids Taxa ids to be extracted, if `NULL`, all taxon in `umi` will be
#' counting.
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @export
taxa_counts <- function(umi, taxids = NULL) {
    if (!is.null(taxids)) {
        umi <- umi$filter(pl$col("taxid")$is_in(taxids))
    }
    # create barcode umi data -----------------------------------
    counts <- umi$lazy()$
        with_columns(pl$col(pl$List(pl$String))$list$join("|"))$
        group_by("barcode", "taxid", "taxa", "rank", "taxon", "ranks")$
        agg(pl$col("umi")$n_unique())$
        with_columns(pl$col("ranks", "taxon")$str$split("|"))$
        #     with_columns(
        #     pl$col("ranks", "taxon")$list$eval(
        #         pl$int_range(0L, pl$element()$len())$add(1L)
        #     )$name$suffix("_len")
        # )$
        #     explode("^.+_len$")$
        #     with_columns(
        #     pl$col("ranks")$list$slice(0L, length = pl$col("ranks_len")),
        #     pl$col("taxon")$list$slice(0L, length = pl$col("taxon_len"))
        # )$
        #     drop("^.+_len$")$
        with_row_index("index")$
        explode(pl$col("ranks", "taxon"))$
        with_columns(
        # A rank code, indicating (U)nclassified, (R)oot, (D)omain,
        # (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or
        # (S)pecies. Taxa that are not at any of these 10 ranks have a rank
        # code that is formed by using the rank code of the closest ancestor
        # rank with a number indicating the distance from that rank. E.g.,
        # "G2" is a rank code indicating a taxon is between genus and
        # species and the grandparent taxon is at the genus rank.
        pl$col("ranks")$
            str$replace("^R", "root")$
            str$replace("^D", "domain")$
            str$replace("^K", "kingdom")$
            str$replace("^P", "phylum")$
            str$replace("^C", "class")$
            str$replace("^O", "order")$
            str$replace("^F", "family")$
            str$replace("^G", "genus")$
            str$replace("^S", "species")
    )$
        collect()$
        pivot(
        values = "taxon",
        index = c("index", "barcode", "taxid", "taxa", "rank", "umi"),
        columns = "ranks"
    )$
        drop("index")
    # relocate columns ----------------------------
    cs <- list(pl$col("barcode", "taxid", "taxa", "rank"))
    columns <- counts$columns
    for (taxa in c("root", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")) {
        pattern <- sprintf("^%s\\d*$", taxa)
        if (any(grepl(pattern, columns, perl = TRUE))) {
            cs <- c(cs, list(pl$col(pattern)))
        }
    }
    counts$select(c(cs, list(pl$col("umi")$alias("counts"))))
}
