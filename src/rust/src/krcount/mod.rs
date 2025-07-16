use anyhow::anyhow;
use anyhow::Result;
use bytes::Bytes;
use extendr_api::prelude::*;
use rustc_hash::FxHashMap as HashMap;
use rustc_hash::FxHashSet as HashSet;

mod count;

use crate::kreport::taxonomy_kreport;
use crate::utils::*;

#[extendr]
fn krcount(
    koutreads: &str,
    kreport: &str,
    umi_tag: Option<&str>,
    barcode_tag: Option<&str>,
    taxonomy: Robj,
    batch_size: usize,
    nqueue: Option<usize>,
) -> std::result::Result<List, String> {
    krcount_internal(
        koutreads,
        kreport,
        umi_tag,
        barcode_tag,
        taxonomy,
        batch_size,
        nqueue,
    )
    .map_err(|e| format!("{}", e))
}

fn krcount_internal(
    koutreads: &str,
    kreport: &str,
    umi_tag: Option<&str>,
    barcode_tag: Option<&str>,
    taxonomy: Robj,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<List> {
    let kreports = taxonomy_kreport(kreport, taxonomy)?;

    // ─── Build taxonomic ancestry map ───────────────────
    // Each taxid maps to a set of its ancestor taxids (inclusive)
    let taxid_to_ancestors = kreports
        .iter()
        .map(|report| {
            // Each report's `taxids` field holds the lineage (ancestors) of this taxon
            let ancestors = report
                .taxids
                .iter()
                .map(|x| x.as_slice())
                .collect::<HashSet<&[u8]>>();
            // Map: current taxid → set of its ancestors
            (report.taxid.as_slice(), ancestors)
        })
        .collect::<HashMap<&[u8], HashSet<&[u8]>>>();

    // ─── Count reads and kmers per (barcode, taxon) ─────
    let counts_map = count::count_kmers_and_reads(
        koutreads,
        taxid_to_ancestors,
        umi_tag,
        barcode_tag,
        batch_size,
        nqueue,
    )?;

    // ─── Determine all observed rank codes ───────────────
    // Examples: U, R, D, K, P, C, O, F, G, S, G2, S1, etc.
    // we first extract all rank codes
    // A rank code, indicating (U)nclassified, (R)oot, (D)omain,
    // (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or
    // (S)pecies. Taxa that are not at any of these 10 ranks have a rank
    // code that is formed by using the rank code of the closest ancestor
    // rank with a number indicating the distance from that rank. E.g.,
    // "G2" is a rank code indicating a taxon is between genus and
    // species and the grandparent taxon is at the genus rank.
    let rank_sets = kreports
        .iter()
        .map(|report| report.ranks.as_slice())
        .flatten()
        .collect::<HashSet<_>>();
    let mut ordered_ranks: Vec<_> = rank_sets.into_iter().collect();
    ordered_ranks.sort_by_key(|r| rank_order_key(r));
    let taxa_cols = ordered_ranks
        .iter()
        .map(|bytes| unsafe { String::from_utf8_unchecked(bytes.to_vec()) })
        .collect::<Vec<_>>();

    // ─── Build taxa table by rank (columns) ──────────────
    // taxa_table: rank → vector of taxon names per report
    let mut taxa_table: HashMap<&[u8], Vec<Rstr>> =
        HashMap::with_capacity_and_hasher(ordered_ranks.len(), rustc_hash::FxBuildHasher);
    for &rank in &ordered_ranks {
        let mut taxa_vec = Vec::with_capacity(kreports.len());
        for report in &kreports {
            let mut value = Rstr::na();
            for (r, t) in report.ranks.iter().zip(report.taxa.iter()) {
                if r == rank {
                    value = u8_to_rstr(t.clone());
                    break;
                }
            }
            taxa_vec.push(value);
        }
        taxa_table.insert(rank, taxa_vec);
    }
    let taxa_vec = ordered_ranks
        .iter()
        .filter_map(|rank| taxa_table.remove(rank.as_slice()))
        .collect::<Vec<_>>();

    // ─── Build data tables: taxon x barcode stats ────────
    // Each table holds rows for barcodes, columns for taxa
    let barcodes = counts_map.keys().into_iter().collect::<Vec<_>>();
    let mut counts_table: HashMap<&Bytes, Vec<Option<usize>>> =
        HashMap::with_capacity_and_hasher(barcodes.len(), rustc_hash::FxBuildHasher);
    let mut kmer_total_table = counts_table.clone();
    let mut kmer_unique_table = counts_table.clone();
    for &barcode in &barcodes {
        let mut reads_vec = Vec::with_capacity(kreports.len());
        let mut kmer_total_vec = Vec::with_capacity(kreports.len());
        let mut kmer_unique_vec = Vec::with_capacity(kreports.len());
        for report in &kreports {
            if let Some(barcode_map) = counts_map.get(barcode) {
                if let Some(reads_and_kmer) = barcode_map.get(report.taxid.as_slice()) {
                    reads_vec.push(Some(reads_and_kmer.reads()));
                    kmer_total_vec.push(Some(reads_and_kmer.kmer_total()));
                    kmer_unique_vec.push(Some(reads_and_kmer.kmer_unique()));
                    continue;
                }
            }
            reads_vec.push(None);
            kmer_total_vec.push(None);
            kmer_unique_vec.push(None);
        }
        counts_table.insert(barcode, reads_vec);
        kmer_total_table.insert(barcode, kmer_total_vec);
        kmer_unique_table.insert(barcode, kmer_unique_vec);
    }
    let counts_vec = barcodes
        .iter()
        .filter_map(|barcode| counts_table.remove(*barcode))
        .collect::<Vec<_>>();
    let kmer_total_vec = barcodes
        .iter()
        .filter_map(|barcode| kmer_total_table.remove(*barcode))
        .collect::<Vec<_>>();
    let kmer_unique_vec = barcodes
        .iter()
        .filter_map(|barcode| kmer_unique_table.remove(*barcode))
        .collect::<Vec<_>>();
    let barcode_cols = barcodes
        .iter()
        .map(|bytes| unsafe { String::from_utf8_unchecked(bytes.to_vec()) })
        .collect::<Vec<_>>();

    Ok(list![
        taxa = List::from_names_and_values(taxa_cols, taxa_vec)
            .map_err(|e| anyhow!("Failed to create list for taxa: {}", e))?,
        counts = List::from_names_and_values(barcode_cols.clone(), counts_vec)
            .map_err(|e| anyhow!("Failed to create list for counts: {}", e))?,
        kmer_total = List::from_names_and_values(barcode_cols.clone(), kmer_total_vec)
            .map_err(|e| anyhow!("Failed to create list for kmer_total: {}", e))?,
        kmer_unique = List::from_names_and_values(barcode_cols, kmer_unique_vec)
            .map_err(|e| anyhow!("Failed to create list for kmer_unique: {}", e))?,
    ])
}

fn rank_order_key(rank: &[u8]) -> (usize, usize) {
    match rank {
        b"U" => (0, 0),
        b"R" => (1, 0),
        b"D" => (2, 0),
        b"K" => (3, 0),
        b"P" => (4, 0),
        b"C" => (5, 0),
        b"O" => (6, 0),
        b"F" => (7, 0),
        b"G" => (8, 0),
        b"S" => (9, 0),
        _ if rank.len() > 1 => (
            match unsafe { rank.get_unchecked(0) } {
                b'U' => 0,
                b'R' => 1,
                b'D' => 2,
                b'K' => 3,
                b'P' => 4,
                b'C' => 5,
                b'O' => 6,
                b'F' => 7,
                b'G' => 8,
                b'S' => 9,
                _ => 10,
            },
            parse_usize(&rank[1 ..]).map_or_else(|_| usize::max_value(), |x| x),
        ),
        _ => (10, 0),
    }
}

extendr_module! {
    mod krcount;
    fn krcount;
}
