#' Tabulating k-mer statistics
#' @description The function count the number of k-mers and unique k-mers
#' assigned to a taxon across barcodes. The cell barcode and unique molecular
#' identifier (UMI) are used to identify unique barcodes and reads. Data is
#' reported for taxa of pre-specified ranks (default genus + species) taking
#' into account all subsequently higher resolution ranks. Reads with any k-mers
#' mapped to the host (e.g. human) are discarded. Reads with >50% of the k-mers
#' map outside the taxon's lineage are also discarded. The output is a table of
#' barcodes, taxonomic IDs, number of k-mers, and number of unique k-mers.
#' @param fa1,fa2 The path to microbiome fasta 1 and 2 file (returned by
#'   `extract_microbiome`). We'll extract barcode sequence from `fa1`.
#' @param microbiome_out The path of microbiome output file (returned by
#'   [extract_microbiome]). if `NULL`, will locate
#'   __sample.microbiome.output.txt__ in `out_dir`.
#' @param ranks Taxa ranks to analyze.
#' @param cb_len Nucleutide length of cell barcodes
#' @param umi_len Nucleutide length of umis
#' @param host Host taxid to exclude
#' @param nsample Max number of barcodes to sample per taxa.
#' @param kmer_len Kraken kmer length.
#' @param min_frac Minimum fraction of kmers directly assigned to taxid to use
#'   read.
#' @inheritParams extract_microbiome
#' @inheritParams run_kraken2
#' @seealso <https://github.com/sjdlabgroup/SAHMI>
#' @export
run_sckmer <- function(fa1, fa2 = NULL, kraken_report = NULL, mpa_report = NULL, microbiome_out = NULL, sample = NULL, out_dir = getwd(), ranks = c("G", "S"), cb_len = 16L, umi_len = 10L, host = 9606L, nsample = 1000L, kmer_len = 35L, min_frac = 0.5, cores = availableCores()) {
    sample <- sample %||% sub("_0*[12]?\\.fa$", "", basename(fa1), perl = TRUE)
    kraken_report <- define_path(kraken_report,
        sample = sample,
        dir = out_dir
    )
    mpa_report <- define_path(mpa_report, sample = sample, dir = out_dir)
    microbiome_out <- define_path(microbiome_out,
        sample = sample,
        dir = out_dir
    )

    # read in fasta data -----------------------------------------------
    reads1 <- ShortRead::readFasta(fa1)
    sequences1 <- ShortRead::sread(reads1)

    barcode <- substr(sequences1, 1L, cb_len)
    umi <- substr(sequences1, cb_len + 1L, cb_len + umi_len)
    # we use id to match data
    id <- gsub("\\s.*", "", ShortRead::id(reads1), perl = TRUE)

    # only keep data in both sequence ----------------------------------
    if (!is.null(fa2)) {
        reads2 <- ShortRead::readFasta(fa2)
        sequences2 <- ShortRead::sread(reads2)
        id2 <- gsub("\\s.*", "", ShortRead::id(reads2), perl = TRUE)
        inter_id <- intersect(id, id2)
        idx <- id %in% inter_id
        sequences1 <- sequences1[idx]
        barcode <- barcode[idx]
        umi <- umi[idx]
        id <- id[idx]
        sequences2 <- sequences2[id2 %in% inter_id]
    } else {
        sequences2 <- NULL
    }

    # prepare kr, mpa and microbiome_out data ---------------------------
    kr <- data.table::fread(kraken_report, header = FALSE, sep = "\t")[-c(1:2)]
    kr[, V8 := str_trim(gsub("[^[:alnum:]]+", " ", V8, perl = TRUE))] # nolint

    mpa <- data.table::fread(mpa_report, header = FALSE, sep = "\t")
    mpa[, taxid := vapply(strsplit(V1, "|", fixed = TRUE), function(x) { # nolint
        str <- sub(
            ".*__", "",
            unlist(x, recursive = FALSE, use.names = FALSE),
            perl = TRUE
        )
        str <- str_trim(gsub("[^[:alnum:]]+", " ", str, perl = TRUE))
        paste0("*", paste0(
            kr$V7[data.table::chmatch(str, kr$V8)],
            collapse = "*"
        ), "*")
    }, character(1L))]
    # mpa$taxid[1L] <- NA_character_
    microbiome_out <- data.table::fread(microbiome_out,
        header = FALSE, drop = 1L
    )
    microbiome_out[
        , c("name", "taxid") := {
            x <- str_match(V3, "\\s*(.+)\\s*\\(taxid\\s*(\\d+)\\s*\\)") # nolint
            apply(x[, -1L, drop = FALSE], 2L, identity, simplify = FALSE)
        }
    ]
    microbiome_out[, taxid := as.integer(taxid)] # nolint

    # extract all necessary taxa -------------------------------------------
    tx <- setdiff(kr$V7[kr$V6 %in% ranks], host) # nolint
    tx <- intersect(tx, microbiome_out$taxid)

    # define kmer for each tx --------------------------------------------
    out <- define_kmer(
        taxa_vec = tx, mpa_report = mpa,
        microbiome_out = microbiome_out,
        host = host, id = id, barcode = barcode, umi = umi,
        sequences1 = sequences1, sequences2 = sequences2,
        nsample = nsample, kmer_len = kmer_len, min_frac = min_frac,
        cores = cores
    )
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    data.table::fwrite(out,
        file = file_path(out_dir, sample, ext = "sckmer.txt"),
        sep = "\t", row.names = FALSE, col.names = TRUE
    )
    invisible(out)
}
utils::globalVariables(c("V8", "taxid", "V1", "V3"))

define_kmer <- function(taxa_vec, mpa_report, microbiome_out, host, id, barcode, umi, sequences1, sequences2, nsample, kmer_len, min_frac, cores) {
    cli::cli_alert("Defining kmer for {.val {length(taxa_vec)}} taxa")
    old_handlers <- new_handlers()
    on.exit(progressr::handlers(old_handlers))
    old_plan <- future::plan("multicore", workers = cores)
    on.exit(future::plan(old_plan), add = TRUE)
    p <- progressr::progressor(along = taxa_vec, auto_finish = TRUE)
    barcode_kmer_list <- future.apply::future_lapply(taxa_vec, function(taxa) {
        # we will exit this function in the median of the process
        # so we update progress in the beginning
        p()
        full_taxa <- grep(paste0("\\*", taxa, "\\*"),
            mpa_report$taxid, # nolint
            perl = TRUE, value = TRUE
        )
        full_taxa <- lapply(strsplit(full_taxa, "*", fixed = TRUE), function(x) {
            x <- x[x != "NA" & x != "" & !is.na(x)]
            as.integer(x)
        })
        child_taxa <- lapply(full_taxa, function(x) {
            x[cumsum(x == taxa) > 0L]
        })
        full_taxa <- unique(unlist(
            full_taxa,
            recursive = FALSE, use.names = FALSE
        ))
        child_taxa <- unique(unlist(
            child_taxa,
            recursive = FALSE, use.names = FALSE
        ))

        # extract output data ------------------------------------------
        taxa_out <- microbiome_out[taxid %in% child_taxa] # nolint
        taxa_out[, c("r1", "r2") := data.table::tstrsplit(V5, # nolint
            split = "|:|", fixed = TRUE
        )]

        taxa_out <- taxa_out[{
            .idx <- sapply(list(r1, r2), function(.x) { # nolint
                !(grepl(paste0(" ", host, ":"), .x, fixed = TRUE) | .x == "")
            })
            if (is.matrix(.idx)) {
                rowSums(.idx, na.rm = TRUE) > 0L
            } else {
                sum(.idx, na.rm = TRUE) > 0L
            }
        }]
        taxa_out[, c("r1", "r2") := lapply(.SD, str_trim), # nolint
            .SDcols = c("r1", "r2")
        ]

        if (nrow(taxa_out) == 0L) {
            return(NULL)
        }

        # extract barcode and umi data ---------------------------------
        idx <- which(id %in% taxa_out$V2)
        taxa_barcode <- barcode[idx]
        taxa_umi <- umi[idx]
        taxa_sequences1 <- as.character(sequences1[idx])
        taxa_sequences2 <- as.character(sequences2[idx])

        if (nrow(taxa_out) > nsample) {
            n <- sample.int(nrow(taxa_out), size = nsample)
        } else {
            n <- seq_len(nrow(taxa_out))
        }
        taxa_out <- taxa_out[n]
        taxa_barcode <- taxa_barcode[n]
        taxa_umi <- taxa_umi[n]
        taxa_sequences1 <- taxa_sequences1[n]
        taxa_sequences2 <- taxa_sequences2[n]

        out_list <- lapply(c("r1", "r2"), function(mate) {
            taxid_nkmer_pair_list <- strsplit(
                taxa_out[[mate]],
                split = "\\s", perl = TRUE
            )
            mate_seq <- switch(mate,
                r1 = taxa_sequences1,
                r2 = taxa_sequences2
            )
            if (all(is.na(mate_seq)) || is.null(mate_seq)) {
                return(NULL)
            }
            barcode_data <- .mapply(function(pair, sequence, barcode) {
                if (length(pair) == 0L) {
                    return(NULL)
                }
                res <- data.table::tstrsplit(
                    pair,
                    split = ":", fixed = TRUE,
                    type.convert = TRUE
                )
                data.table::setDT(res)
                data.table::setnames(res, c("taxid", "nkmer"))
                res[, c("fkmer", "nt_start", "nt_end") := list(
                    nkmer / sum(nkmer, na.rm = TRUE), # nolint
                    cumsum(nkmer) - nkmer + 1L,
                    cumsum(nkmer) + kmer_len - 1L
                )]
                res[, nt_len := nt_end - nt_start + 1L] # nolint
                res <- res[taxid %in% c(0L, full_taxa)] # nolint
                if (nrow(res) == 0L || is.na(sequence) || sum(res$fkmer, na.rm = TRUE) < min_frac) {
                    return(NULL)
                }
                # calculate kmer
                kmer <- .mapply(function(k, start) {
                    m <- seq_len(k)
                    substring(
                        sequence,
                        start + m - 1L,
                        start + m + kmer_len - 2L
                    )
                }, unname(res[, c("nkmer", "nt_start")]), NULL)
                data.table::data.table(
                    barcode = barcode,
                    taxid = taxa,
                    k = unlist(kmer, recursive = FALSE, use.names = FALSE),
                    n = sum(res$nt_len[res$taxid %in% child_taxa], na.rm = TRUE)
                )
            }, list(
                pair = taxid_nkmer_pair_list,
                sequence = mate_seq,
                barcode = taxa_barcode
            ), NULL)
            data.table::rbindlist(barcode_data, use.names = FALSE)
        })
        out_data <- data.table::rbindlist(out_list, use.names = FALSE)
        if (nrow(out_data) == 0L) {
            return(NULL)
        }
        out_data[, list(
            kmer = length(k), uniq = length(unique(k)) # nolint
        ), by = c("barcode", "taxid")]
    })
    data.table::rbindlist(barcode_kmer_list, use.names = FALSE)
}
utils::globalVariables(
    c("taxid", "V5", "r1", "r2", "nkmer", "nt_len", "nt_end", "nt_start", "k")
)
