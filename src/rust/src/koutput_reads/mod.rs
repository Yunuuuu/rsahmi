use std::fs::File;

use aho_corasick::{AhoCorasick, AhoCorasickKind};
use anyhow::{anyhow, Result};
use rustc_hash::FxHashMap as HashMap;
use rustc_hash::FxHashSet as HashSet;

mod koutput;
// mod reads;

use crate::kreport::parse_kreport;
use crate::reader::bytes::BytesProgressBarReader;

fn koutput_reads(
    kreport: &str,
    koutput: &str,
    fq: &str,
    taxonomy: Option<Vec<&str>>,
    lca: Option<Vec<&str>>,
    exclude: Option<Vec<&str>>,
    polyn_threshold: usize,
    phred_threshould: usize,
    ofile: &str,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<()> {
    let mut kreports = parse_kreport(kreport)?;

    if let Some(taxonomy) = taxonomy {
        // Parse taxon strings like "rank__name" into rank-name pairs
        let rank_taxon_sets = taxonomy
            .iter()
            .filter_map(|t| {
                let mut pair = t.splitn(2, "__");
                if let (Some(rank), Some(taxa)) = (pair.next(), pair.next()) {
                    Some((rank.as_bytes(), taxa.as_bytes()))
                } else {
                    None
                }
            })
            .collect::<HashSet<(&[u8], &[u8])>>();

        // Fail early if no valid taxon entries
        if !taxonomy.is_empty() && rank_taxon_sets.is_empty() {
            return Err(anyhow!("No valid taxonomy provided. 'taxonomy' must be in the format 'rank__name', where 'rank' and 'name' are separated by '__'."));
        }

        // Parsing kraken2 report: only contain information specified by `taxon`
        kreports = kreports
            .into_iter()
            .filter(|kr| {
                kr.ranks
                    .iter()
                    .zip(kr.taxa.iter())
                    .any(|(rank, taxa)| rank_taxon_sets.contains(&(rank, taxa)))
            })
            .collect();
    }

    // Build a map: taxid → set of its ancestor taxids
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

    // For each taxid, find all other taxids that consider it an ancestor
    // i.e., build a reverse lookup: taxid → set of descendant taxids
    let taxid_to_descendants = kreports
        .iter()
        .map(|report| {
            let taxid = report.taxid.as_slice();
            // Iterate over all taxids and check whose ancestor set includes the current taxid
            let descendants = taxid_to_ancestors
                .iter()
                .filter_map(|(child, parents)| {
                    if parents.contains(taxid) {
                        Some(*child)
                    } else {
                        None
                    }
                })
                .collect::<HashSet<&[u8]>>();

            // Map: current taxid → set of its descendant taxids
            (taxid, descendants)
        })
        .collect::<HashMap<&[u8], HashSet<&[u8]>>>();

    let patterns = kreports
        .iter()
        // Always include the descendants
        .filter_map(|kr| taxid_to_descendants.get(kr.taxid.as_slice()))
        .flatten()
        .map(|taxid| {
            let mut v = Vec::with_capacity(b"(taxid ".len() + taxid.len() + 1); // estimated capacity
            v.extend_from_slice(b"(taxid ");
            v.extend_from_slice(taxid);
            v.push(b')');
            v
        })
        .collect::<Vec<_>>();

    let include_aho = AhoCorasick::builder()
        .kind(Some(AhoCorasickKind::DFA))
        .build(patterns)?;

    // A space-delimited list indicating the LCA mapping of each
    // k-mer in the sequence(s). For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
    //
    // the first 13 k-mers mapped to taxonomy ID #562
    // the next 4 k-mers mapped to taxonomy ID #561
    // the next 31 k-mers contained an ambiguous nucleotide
    // the next k-mer was not in the database
    // the last 3 k-mers mapped to taxonomy ID #562
    let exclude_aho = exclude
        .map(|v| {
            let patterns: Vec<Vec<u8>> = v
                .iter()
                .map(|taxid| {
                    let mut pattern = Vec::with_capacity(taxid.len() + 1); // estimated capacity
                    pattern.extend_from_slice(taxid.as_bytes());
                    pattern.push(b':');
                    pattern
                })
                .collect();
            AhoCorasick::builder()
                .kind(Some(AhoCorasickKind::DFA))
                .build(patterns)
        })
        .transpose()?;

    // Read Kraken2 output and extract matched records
    let koutput_file = File::open(koutput)?;
    let reader = BytesProgressBarReader::new(koutput_file);
    let koutput = koutput::reader_parse_koutput(
        reader,
        include_aho,
        exclude_aho,
        chunk_size,
        batch_size,
        nqueue,
    )?;

    // For each koutput row, we calculate kmer information
    Ok(())
}
