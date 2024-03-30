#' Tabulating k-mer statistics
#'
#' @description The function count the number of k-mers and unique k-mers
#' assigned to a taxon across barcodes. The cell barcode and unique molecular
#' identifier (UMI) are used to identify unique barcodes and reads. Data is
#' reported for taxa of pre-specified ranks (default genus + species) taking
#' into account all subsequently higher resolution ranks. The output is a table
#' of barcodes, taxonomic IDs, number of k-mers, and number of unique k-mers.
#'
#' @param fa1,fa2 The path to microbiome fasta 1 and 2 file (returned by
#'   [extract_kraken_reads]).
#' @inheritParams extractor
#' @param kraken_out The path of microbiome output file. Usually should be
#'   filtered with [extract_kraken_output].
#' @param cb_and_umi A function takes sequence id, read1, read2 and return a
#' list of 2 corresponding to cell barcode and UMI respectively., each should
#' have the same length of the input.
#' @param ranks Taxa ranks to analyze.
#' @param exclude A character of taxid to exclude, for `SAHMI`, the host taxid.
#' Reads with any k-mers mapped to the `exclude` are discarded.
#' @param kmer_len Kraken kmer length. Default: `35L`, which is the default kmer
#' size of kraken2.
#' @param min_frac Minimum fraction of kmers directly assigned to taxid to use
#' read. Reads with `<=min_frac` of the k-mers map inside the taxon's lineage
#' are also discarded.
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @examples
#' # for sequence from umi-tools, we can use following function
#' cb_and_umi <- function(sequence_id, read1, read2) {
#'     out <- lapply(
#'         strsplit(sequence_id, "_", fixed = TRUE),
#'         `[`, 2:3
#'     )
#'     lapply(1:2, function(i) {
#'         vapply(out, function(o) as.character(.subset2(o, i)), character(1L))
#'     })
#' }
#' @export
#' @importFrom polars pl
sckmer <- function(fa1, kraken_report, kraken_out, fa2 = NULL,
                   cb_and_umi = function(sequence_id, read1, read2) {
                       list(
                           substring(read1, 1L, 16L),
                           substring(read1, 17L, 28L)
                       )
                   },
                   ranks = c("G", "S"), kmer_len = 35L,
                   min_frac = 0.5, exclude = "9606",
                   threads = 10L) {
    kreport <- parse_kraken_report(kraken_report)
    exclude <- pl$Series(values = exclude)$cast(pl$String)

    # extract operated taxon -----------------------------------------
    taxon_struct <- kreport$filter(pl$col("ranks")$list$last()$is_in(ranks))$
        select(
        pl$col("taxids")$list$last()$alias("taxid"),
        pl$col("taxon")$list$last()$alias("taxa")
    )$
        filter(pl$col("taxid")$is_in(exclude)$not())$
        to_struct()

    # we get all operated taxon and their children
    # this is just used to filter kraken output data
    # otherwise, kraken output data will be very large
    children_taxon <- taxon_children(
        kreport,
        taxon_struct$struct$field("taxid")
    )

    # prepare taxid:kmer data ------------------------------------------
    kout <- pl$scan_csv(kraken_out, has_header = FALSE, separator = "\t")$
        filter(
        pl$col("column_5")$str$
            contains_any(pl$concat_str(pl$Series(values = ":"), exclude))$
            not()
    )$
        select(
        pl$col("column_2")$alias("sequence_id"),
        pl$col("column_3")$str$
            extract(pl$lit("\\s*(.+?)\\s*\\(taxid\\s*(\\d+|A)\\s*\\)"), 1L)$
            alias("name"),
        pl$col("column_3")$str$
            extract(pl$lit("\\s*(.+?)\\s*\\(taxid\\s*(\\d+|A)\\s*\\)"), 2L)$
            alias("taxid"),
        # Note that paired read data will contain a "|:|" token in this list to
        # indicate the end of one read and the beginning of another.
        pl$col("column_5")$str$split("|:|")$alias("LCA")
    )$
        filter(pl$col("taxid")$is_in(children_taxon))$
        with_row_index("index")$
        explode("LCA")$ # all reads have been included in one column
        with_columns(
        # split LCA column into taxid:kmer pairs
        # A space-delimited list indicating the LCA mapping of each k-mer in the
        # sequence(s). For example, "562:13 561:4 A:31 0:1 562:3" would indicate
        # that:
        # the first 13: k-mers mapped to taxonomy ID #562
        # the next 4: k-mers mapped to taxonomy ID #561
        # the next 31: k-mers contained an ambiguous nucleotide
        # the next k-mer was not in the database
        # the last 3 k-mers mapped to taxonomy ID #562
        pl$col("LCA")$str$strip_chars(),
        pl$concat_str(
            pl$lit("read"),
            pl$int_range(end = pl$len())$over("index")$add(1L)$cast(pl$String),
            separator = ""
        )$alias("header")
    )$
        collect()$
        # pivot method must run with DataFrame
        pivot("LCA",
        index = c("index", "name", "taxid", "sequence_id"),
        columns = "header"
    )

    # check read ----------------------------------------------------
    read_nms <- setdiff(
        kout$columns,
        c("index", "name", "taxid", "sequence_id")
    )
    if (length(read_nms) == 1L) {
        if (!is.null(fa2)) {
            cli::cli_warn(paste(
                "{.arg fa2} will be ignored",
                "only one read was found in {.field kraken}",
                sep = ", "
            ))
            fa2 <- NULL
        }
    } else if (length(read_nms) == 2L) {
        if (is.null(fa2)) {
            cli::cli_warn(paste(
                "read2 in {.arg kraken_report} will be ignored",
                "since {.arg fa2} was not provided"
            ))
            read_nms <- "read1"
            kout <- kout$select(
                pl$col("index", "name", "taxid", "sequence_id", read_nms)
            )
        }
    } else {
        cli::cli_abort("Invalid kraken2 output file")
    }

    # extract kmer information  -------------------------------------
    kout <- kout$lazy()$with_columns(
        pl$col(read_nms)$str$extract_all("(\\d+|A):")$name$suffix("_taxid"),
        pl$col(read_nms)$str$extract_all(":\\d+")$name$suffix("_kmer")
    )$
        explode(pl$col("^.+_taxid$", "^.+_kmer$"))$
        with_columns(
        pl$col("^.+_taxid$")$str$strip_chars_end(":"),
        pl$col("^.+_kmer$")$str$strip_chars_start(":")$cast(pl$Int64)
    )$
        group_by("index", "name", "taxid", "sequence_id", read_nms,
        maintain_order = FALSE
    )$
        agg(pl$col("^.+_kmer$", "^.+_taxid$"))$
        with_columns(
        pl$col("^.+_kmer$")$list$
            eval(pl$element()$div(pl$element()$sum()$cast(pl$Float64)))$
            name$suffix("_frac"),
        pl$col("^.+_kmer$")$list$
            eval(pl$element()$cum_sum()$sub(pl$element())$add(1L))$
            name$suffix("_nt_start"),
        pl$col("^.+_kmer$")$list$
            eval(pl$element()$cum_sum()$add(kmer_len)$sub(1L))$
            name$suffix("_nt_end"),
        pl$col("^.+_kmer$")$list$
            eval(pl$element()$add(kmer_len)$sub(1L))$
            name$suffix("_nt_len")
    )$
        collect()

    # read in fasta data -----------------------------------------------
    read1 <- ShortRead::readFasta(fa1)
    id1 <- as.character(ShortRead::id(read1))
    read1 <- ShortRead::sread(read1)

    if (!is.null(fa2)) {
        read2 <- ShortRead::readFasta(fa2)
        id2 <- as.character(ShortRead::id(read2))
        read2 <- ShortRead::sread(read2)

        # only keep data in both sequence ----------------------------------
        ids <- intersect(id1, id2)
    } else {
        ids <- id1
        read2 <- NULL
    }

    # only keep sequence in fa1, fa2, and kraken output
    kout <- kout$filter(pl$col("sequence_id")$is_in(ids))
    ids <- kout$get_column("sequence_id")$to_r()
    read1 <- read1[match(ids, id1)]
    if (!is.null(fa2)) read2 <- read2[match(ids, id2)]

    # for operated taxon, we also remove items not in kraken output
    taxon_struct <- taxon_struct$
        filter(
        taxon_struct$struct$field("taxid")$is_in(kout$get_column("taxid"))
    )$
        unique()

    # extract cell barcode and umi -----------------------------------
    cb_and_umi <- cb_and_umi(ids, read1, read2)
    if (length(cb_and_umi) != 2L) {
        cli::cli_abort(c(
            "{.code length(cb_and_umi) == 2L} is not {.code TRUE}",
            i = "{.fn cb_and_umi} must return a list of length 2"
        ))
    } else if (!all(lengths(cb_and_umi) == length(ids))) {
        cli::cli_abort(c(
            paste(
                "{.code all(lengths(cb_and_umi) == length(sequence_id))}",
                "is not {.code TRUE}"
            ),
            i = paste(
                "{.fn cb_and_umi} must return cell barcode and umi",
                "for each reads"
            )
        ))
    }

    # integrate sequence, cell barcode and umi ----------------------
    kout <- kout$with_columns(
        read1_sequence = pl$Series(values = as.character(read1))
    )
    if (!is.null(fa2)) {
        kout <- kout$with_columns(
            read2_sequence = pl$Series(values = as.character(read2))
        )
    }

    # every row correspond to a single sequence read
    kout <- kout$with_columns(
        cb = pl$Series(values = cb_and_umi[[1L]]),
        umi = pl$Series(values = cb_and_umi[[2L]])
    )

    # define kmer ---------------------------------------------------
    kmer_list <- polars_lapply(
        taxon_struct, kmer_query,
        kout = kout, kreport = kreport,
        read_nms = read_nms, kmer_len = kmer_len,
        min_frac = min_frac,
        .progress = list(
            name = "Defining kmer",
            format = paste(
                "{cli::pb_bar} {cli::pb_current}",
                "{cli::pb_total} [{cli::pb_rate}] | {cli::pb_eta_str}",
                sep = "/"
            ),
            format_done = paste(
                "{cli::pb_name} for {.val {cli::pb_total}} tax{?a/on}",
                "in {cli::pb_elapsed}"
            )
        ),
        .threads = threads
    )
    kmer <- pl$concat(kmer_list, how = "vertical")$
        select(
        pl$col("cb")$alias("barcode"),
        pl$col("taxid", "taxa", "kmer_len", "kmer_n_unique")
    )

    # prepare data for taxa counting by UMI ----------------------
    umi <- kreport$select(
        pl$col("taxids")$list$last()$alias("taxid"),
        pl$col("taxon")$list$last()$alias("taxa"),
        pl$col("ranks")$list$last()$alias("rank"),
        pl$col("taxon"), pl$col("ranks")
    )$join(
        kout$select(
            pl$col("cb")$alias("barcode"),
            pl$col("taxid", "umi")
        )$unique(),
        on = "taxid", how = "inner"
    )
    list(kmer = kmer, umi = umi)
}

taxon_children <- function(kreport, taxon) {
    kreport$select(
        pl$col("taxids")$list$gather(
            pl$col("taxids")$list$eval(
                pl$arg_where(pl$element()$is_in(taxon)$cum_sum()$gt(0L))
            )
        )$explode()
    )$
        filter(pl$col("taxids")$is_not_null())$
        to_series()$unique()
}

kmer_query <- function(kout, kreport, read_nms, taxa_struct, kmer_len, min_frac) {
    taxa <- taxa_struct$struct$field("taxa")
    taxid <- taxa_struct$struct$field("taxid")
    lineage_report <- kreport$filter(pl$col("taxids")$list$contains(taxid))
    child_taxon <- taxon_children(lineage_report, taxid)
    lineage_taxon <- lineage_report$
        get_column("taxids")$explode()$append("0")$unique()

    lazy_kout <- kout$lazy()$filter(pl$col("taxid")$is_in(child_taxon))
    out_list <- lapply(read_nms, function(read_nm) {
        cols <- c("cb", paste0(read_nm, "_sequence"))
        lazy_kout$
            select(
            pl$col(cols),
            pl$col(paste0("^", read_nm, "_kmer.*$"))$list$gather(
                pl$col(paste0(read_nm, "_taxid"))$list$eval(
                    pl$arg_where(pl$element()$is_in(lineage_taxon))
                )
            )
        )$
            filter(
            pl$col(paste0(read_nm, "_kmer_frac"))$list$sum()$gt_eq(min_frac)
        )$
            # https://github.com/sjdlabgroup/SAHMI/blob/main/functions/sckmer.r#L161
            # official SAHMI define following field but don't use it
            #     with_columns(
            #     pl$col(paste0(read_nm, "_kmer_nt_len"))$list$gather(
            #         pl$col(paste0(read_nm, "_taxid"))$list$eval(
            #             pl$arg_where(pl$element()$is_in(child_taxon))
            #         )
            #     )$list$sum()$alias(paste0(read_nm, "_n"))
            # )$
            select(pl$col(cols, "^read\\d+_kmer(_nt_start)?$"))$
            explode("^read\\d+_kmer(_nt_start)?$")$
            with_columns(
            pl$int_ranges(0L, pl$col("^read\\d_kmer$"))$alias("start")
        )$
            explode("start")$
            select(
            pl$col("cb"),
            pl$col(paste0(read_nm, "_sequence"))$
                str$slice(
                pl$col("start")$add(pl$col(paste0(read_nm, "_kmer_nt_start"))),
                kmer_len
            )$alias("kmer")
        )
    })
    pl$concat(out_list, how = "vertical")$
        # for a sample, every cell barcode is unique, we can just group_by `cb`
        # which is what official SAHMI do.
        group_by("cb", maintain_order = FALSE)$
        agg(
        pl$col("kmer")$len()$alias("kmer_len"),
        pl$col("kmer")$n_unique()$alias("kmer_n_unique")
    )$
        with_columns(taxid = taxid, taxa = taxa)$
        collect(collect_in_background = TRUE)
}
