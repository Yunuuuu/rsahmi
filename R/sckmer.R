#' Quantify K-mers and UMI Data by Taxon
#'
#' This function counts total and unique k-mers per taxon across cell barcodes,
#' using both the cell barcode and unique molecular identifier (UMI) to resolve
#' read identity at the single-cell level. It aggregates k-mer counts for each
#' taxonomic rank of interest (by default, genus and species), including all
#' descendant taxa within those ranks.
#'
#' @inheritParams kractor
#' @param koutput Path to the extracted Kraken2 output file, typically
#'   filtered using [`kractor()`].
#' @param reads A character vector of extracted FASTA files (from
#' [`kractor()`]). Can be one file (single-end) or two files (paired-end).
#' @param extract_kreport Path to the file where the extracted Kraken2 report,
#'   matched to `koutput`, will be saved. Defaults to
#'   `"kraken_microbiome_report.txt"`.
#' @param extract_kmer Path to the file where the quantified k-mer data,
#'   derived from `koutput`, will be saved. Defaults to
#'   `"kraken_microbiome_kmer.txt"`.
#' @param extract_umi Path to the file where the quantified UMI data,
#'   derived from `koutput`, will be saved. Defaults to
#'   `"kraken_microbiome_umi.txt"`.
#' @param barcode_extractor,umi_extractor A function takes sequence IDs, read1,
#' read2 and return a character corresponding to cell barcode and UMI
#' respectively, each should have the same length of the input.
#' @param ranks Taxa ranks to analyze.
#' @param kmer_len Kraken kmer length. Default: `35L`, which is the default kmer
#' size of kraken2.
#' @param min_frac Minimum fraction of kmers directly assigned to taxid to use
#' read. Reads with `<=min_frac` of the k-mers map inside the taxon's lineage
#' are also discarded.
#' @param exclude A character of taxid to exclude, for `SAHMI`, the host taxid.
#' Reads with any k-mers mapped to the `exclude` are discarded.
#' @param threads Integer. Number of threads to use. Default to all
#' available threads.
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @examples
#' \dontrun{
#' # 1. `kreport` should be the output of `blit::kraken2()`.
#' # 2. `koutput` should be the extracted kraken2 output by `kractor()`
#' # 3. `reads` should be the output of `kractor()`
#' sckmer(
#'     kreport = "kraken_report.txt",
#'     koutput = "kraken_microbiome_output.txt",
#'     reads = "kraken_microbiome_reads.fa",
#'     # if you use paired reads.
#'     # reads = c(
#'     #   "kraken_microbiome_reads_1.fa",
#'     #   "kraken_microbiome_reads_2.fa"
#'     # ),
#'     odir = NULL
#' )
#' }
#' @return None. This function generates the following files:
#'  - `extract_kreport`: A filtered version of the Kraken2 taxonomic report,
#'   containing only taxa that meet the `ranks` criteria and are observed in
#'   `koutput`. Used by [`rpmm_quantile()`] and [`slsd()`].
#'  - `extract_kmer`: A table quantifying total and unique k-mers assigned to
#'   each taxon across barcodes. Used by [`blsd()`].
#'  - `extract_umi`: A table of taxon–barcode–UMI combinations indicating all
#'   observed UMI-tagged reads per taxon. Used by [`taxa_counts()`].
#' @export
sckmer <- function(kreport, koutput, reads,
                   extract_kreport = NULL,
                   extract_kmer = NULL,
                   extract_umi = NULL,
                   barcode_extractor = NULL,
                   umi_extractor = NULL,
                   ranks = c("G", "S"), kmer_len = 35L,
                   min_frac = 0.5, exclude = "9606",
                   threads = NULL, odir = NULL) {
    use_polars()
    assert_string(kreport, allow_empty = FALSE)
    assert_string(koutput, allow_empty = FALSE)
    reads <- as.character(reads)
    if (length(reads) == 1L) {
        fa1 <- reads
        fa2 <- NULL
    } else if (length(reads) == 2L) {
        fa1 <- reads[1L]
        fa2 <- reads[2L]
    } else {
        cli::cli_abort("{.arg reads} must be of length 1 or 2")
    }
    assert_string(extract_kreport, allow_empty = FALSE, allow_null = TRUE)
    assert_string(extract_kmer, allow_empty = FALSE, allow_null = TRUE)
    assert_string(extract_umi, allow_empty = FALSE, allow_null = TRUE)
    assert_number_whole(threads,
        min = 1, max = parallel::detectCores(),
        allow_null = TRUE
    )
    threads <- threads %||% parallel::detectCores()
    assert_string(odir, allow_empty = FALSE, allow_null = TRUE)
    if (is.null(odir)) odir <- getwd()
    extract_kreport <- extract_kreport %||% "kraken_microbiome_report.txt"
    extract_kmer <- extract_kmer %||% "kraken_microbiome_kmer.txt"
    extract_umi <- extract_umi %||% "kraken_microbiome_umi.txt"

    cli::cli_alert_info("Parsing {.path {kreport}}")
    kreport <- parse_kraken_report(kreport)
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
    cli::cli_alert_info("Parsing {.path {koutput}}")
    kout <- pl$scan_csv(koutput, has_header = FALSE, separator = "\t")$
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
        # pivot method must run with DataFrame instead of LazyDataFrame
        collect()$
        pivot(
        on = "header",
        values = "LCA",
        index = c("index", "name", "taxid", "sequence_id"),
    )

    # check read ----------------------------------------------------
    read_nms <- setdiff(
        kout$columns,
        c("index", "name", "taxid", "sequence_id")
    )
    if (length(read_nms) == 1L) {
        if (!is.null(fa2)) {
            cli::cli_warn(paste(
                "the second {.arg reads} will be ignored",
                "only one read was found in {.arg koutput}",
                sep = ", "
            ))
            fa2 <- NULL
        }
    } else if (length(read_nms) == 2L) {
        if (is.null(fa2)) {
            cli::cli_warn(paste(
                "read2 in {.arg kreport} will be ignored",
                "since only one {.arg reads} was provided"
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
    )
    for (read_nm in read_nms) {
        kout <- kout$
            explode(pl$col(sprintf(c("%s_taxid", "%s_kmer"), read_nm)))$
            with_columns(
            pl$col(sprintf("%s_taxid", read_nm))$
                str$strip_chars_end(":"),
            pl$col(sprintf("%s_kmer", read_nm))$
                str$strip_chars_start(":")$cast(pl$Int64)
        )$
            group_by(
            "index", "name", "taxid", "sequence_id",
            sprintf(c("%s_taxid", "%s_kmer"), setdiff(read_nms, read_nm)),
            maintain_order = FALSE
        )$
            agg(pl$col(sprintf(c("%s_taxid", "%s_kmer"), read_nm)))
    }
    kout <- kout$
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
    cli::cli_alert_info("Reading the first read: {.path {fa1}}")
    read1 <- ShortRead::readFasta(fa1)
    id1 <- as.character(ShortRead::id(read1))
    id1 <- vapply(strsplit(id1, " |\t"), .subset2, character(1L), 1L)
    read1 <- ShortRead::sread(read1)

    if (!is.null(fa2)) {
        cli::cli_alert_info("Reading the second read: {.path {fa2}}")
        read2 <- ShortRead::readFasta(fa2)
        id2 <- as.character(ShortRead::id(read2))
        id2 <- vapply(strsplit(id2, " |\t"), .subset2, character(1L), 1L)
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
    cli::cli_alert_info("Parsing {.field cell barcode}")
    if (is.null(barcode_extractor)) {
        if (grepl("@RSAHMI:UMI:[a-zA-Z]+:BARCODE:[a-zA-Z]+:RSAHMI@", ids[1L])) {
            barcode <- sub(
                "@RSAHMI:UMI:[a-zA-Z]+:BARCODE:([a-zA-Z]+):RSAHMI@",
                "\\1", ids
            )
        } else {
            cli::cli_abort(c(
                "x" = "No recognizable barcode pattern found in sequence id",
                "!" = "try to provide {.arg barcode_extractor} instead."
            ))
        }
    } else {
        barcode <- as.character(barcode_extractor(ids, read1, read2))
        if (length(barcode) != length(ids)) {
            cli::cli_abort(
                "{.fn barcode_extractor} must return cell barcode for each read"
            )
        }
    }

    cli::cli_alert_info("Parsing {.field UMI}")
    if (is.null(umi_extractor)) {
        if (grepl("@RSAHMI:UMI:[a-zA-Z]+:BARCODE:[a-zA-Z]+:RSAHMI@", ids[1L])) {
            umi <- sub(
                "@RSAHMI:UMI:([a-zA-Z]+):BARCODE:([a-zA-Z]+):RSAHMI@",
                "\\1", ids
            )
        } else {
            cli::cli_abort(c(
                "x" = "No recognizable UMI pattern found in sequence id",
                "!" = "try to provide {.arg umi_extractor} instead."
            ))
        }
    } else {
        umi <- as.character(umi_extractor(ids, read1, read2))
        if (length(umi) != length(ids)) {
            cli::cli_abort("{.fn umi_extractor} must return umi for each read")
        }
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
        cb = pl$Series(values = barcode),
        umi = pl$Series(values = umi)
    )

    # prepare data for blsa ----------------------
    # define kmer ---------------------------------------------------
    cli::cli_alert_info("Calculating {.field kmer}")
    with(
        mirai::daemons(threads, dispatcher = FALSE),
        {
            query_mirai <- mirai::mirai_map( # nolint
                seq_len(taxon_struct$len()),
                kmer_query,
                kout = kout,
                kreport = kreport,
                taxon_struct = taxon_struct,
                read_nms = read_nms,
                kmer_len = kmer_len,
                min_frac = min_frac
            )
            kmer_list <- query_mirai[.progress]
        }
    )
    kmer <- pl$concat(kmer_list, how = "vertical")$
        select(
        pl$col("cb")$alias("barcode"),
        pl$col("taxid", "taxa", "kmer_len", "kmer_n_unique")
    )

    # filter kreport only in kmer data -----------
    kreport <- kreport$
        with_columns(
        pl$col("taxids")$list$last()$alias("taxid"),
        pl$col("taxon")$list$last()$alias("taxa"),
        pl$col("ranks")$list$last()$alias("rank")
    )$
        filter(pl$col("taxid")$is_in(taxon_struct$struct$field("taxid")))

    # prepare data for taxa counting by UMI ----------------------
    # should we include all children taxa for a taxa?
    # SAHMI don't use children taxon
    umi <- kreport$select(
        pl$col("taxid"), pl$col("taxa"), pl$col("rank"),
        pl$col("taxon"), pl$col("ranks")
    )$join(
        kout$select(
            pl$col("cb")$alias("barcode"),
            pl$col("taxid", "umi")
        )$unique(),
        on = "taxid", how = "inner"
    )

    # prepare data for slsd ----------------------
    kreport <- kreport$select(
        pl$col("taxid"), pl$col("taxa"), pl$col("rank"),
        pl$all()$
            exclude(c("taxid", "taxa", "rank", "taxids", "taxon", "ranks"))$
            list$last()
    )

    # save all results ---------------------------
    kreport$write_parquet(file.path(odir, extract_kreport))
    kmer$write_parquet(file.path(odir, extract_kmer))
    umi$write_parquet(file.path(odir, extract_umi))

    cli::cli_inform(c("v" = "Finished"))
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

kmer_query <- function(i, kout, kreport, taxon_struct, read_nms,
                       kmer_len, min_frac) {
    taxa_struct <- taxon_struct$slice(i - 1L, length = 1L)
    taxa <- taxa_struct$struct$field("taxa")
    taxid <- taxa_struct$struct$field("taxid")
    lineage_report <- kreport$filter(pl$col("taxids")$list$contains(taxid))
    child_taxon <- taxon_children(lineage_report, taxid)
    lineage_taxon <- lineage_report$
        get_column("taxids")$explode()$append("0")$unique()

    lazy_kout <- kout$lazy()$filter(pl$col("taxid")$is_in(child_taxon))
    # Note: the SAHMI don't deduplicate UMI, so for each barcode,
    # there may exists duplicates reads?
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
        group_by("cb", maintain_order = FALSE)$
        agg(
        pl$col("kmer")$len()$alias("kmer_len"),
        pl$col("kmer")$n_unique()$alias("kmer_n_unique")
    )$
        with_columns(taxid = taxid, taxa = taxa)$
        collect()
}
