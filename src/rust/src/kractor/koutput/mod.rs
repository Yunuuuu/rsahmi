use aho_corasick::{AhoCorasick, AhoCorasickKind};
use anyhow::{anyhow, Context, Result};
use extendr_api::prelude::*;
use indicatif::{MultiProgress, ProgressBar, ProgressFinish};
use rustc_hash::FxHashMap as HashMap;
use rustc_hash::FxHashSet as HashSet;

use crate::utils::*;

mod parse;

pub(crate) fn kractor_koutput(
    kreport: &str,
    koutput: &str,
    ofile: &str,
    taxonomy: Robj,
    ranks: Robj,
    taxa: Robj,
    taxids: Robj,
    exclude: Robj,
    descendants: bool,
    compression_level: i32,
    batch_size: usize,
    chunk_bytes: usize,
    nqueue: Option<usize>,
    threads: usize,
) -> Result<()> {
    let mut kreports = crate::kreport::parse_kreport(kreport)?;

    let taxonomy =
        robj_to_option_str(&taxonomy).with_context(|| format!("Failed to parse 'taxonomy'"))?;
    let ranks = robj_to_option_str(&ranks).with_context(|| format!("Failed to parse 'ranks'"))?;
    let taxa = robj_to_option_str(&taxa).with_context(|| format!("Failed to parse 'taxa'"))?;
    let taxids =
        robj_to_option_str(&taxids).with_context(|| format!("Failed to parse 'taxids'"))?;
    let exclude =
        robj_to_option_str(&exclude).with_context(|| format!("Failed to parse 'exclude'"))?;

    if taxonomy.is_none()
        && ranks.is_none()
        && taxa.is_none()
        && taxids.is_none()
        && exclude.is_none()
    {
        return Err(anyhow!(
            "One of 'taxonomy', 'ranks', 'taxa', 'taxids', 'exclude' must be provided"
        ));
    }

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
        if kreports.is_empty() {
            return Err(anyhow!(
                "No taxonomic matches found in the kreport file for {:?}.",
                taxonomy
            ));
        }
    }

    let mut targeted_taxids: Vec<&[u8]>;
    if ranks.is_some() || taxa.is_some() || taxids.is_some() {
        // Parse set of desired taxonomic ranks
        let mut reports = kreports.iter().collect::<Vec<_>>();
        if let Some(ranks) = ranks {
            let ranks_sets = ranks
                .iter()
                .map(|x| x.as_bytes())
                .collect::<HashSet<&[u8]>>();
            reports = reports
                .into_iter()
                .filter(|kr| ranks_sets.contains(kr.rank.as_slice()))
                .collect();
        }
        if let Some(taxa) = taxa {
            let taxa_sets = taxa
                .iter()
                .map(|x| x.as_bytes())
                .collect::<HashSet<&[u8]>>();
            reports = reports
                .into_iter()
                .filter(|kr| taxa_sets.contains(kr.taxon.as_slice()))
                .collect();
        }
        if let Some(taxids) = taxids {
            let taxids_sets = taxids
                .iter()
                .map(|x| x.as_bytes())
                .collect::<HashSet<&[u8]>>();
            reports = reports
                .into_iter()
                .filter(|kr| taxids_sets.contains(kr.taxid.as_slice()))
                .collect();
        }
        targeted_taxids = reports.into_iter().map(|kr| kr.taxid.as_slice()).collect();
    } else {
        targeted_taxids = kreports.iter().map(|kr| kr.taxid.as_slice()).collect();
    }

    if descendants {
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
                (report.taxid.as_slice(), descendants)
            })
            .collect::<HashMap<&[u8], HashSet<&[u8]>>>();
        targeted_taxids = targeted_taxids
            .into_iter()
            .filter_map(|taxid| taxid_to_descendants.get(taxid))
            .flatten()
            .copied()
            .collect()
    }

    let include_sets = targeted_taxids.into_iter().collect::<HashSet<&[u8]>>();

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
    let reader_style = progress_reader_style()?;
    let writer_style = progress_writer_style()?;
    let progress = MultiProgress::new();
    let pb1 = progress.add(
        ProgressBar::new(std::fs::metadata(koutput)?.len() as u64)
            .with_finish(ProgressFinish::Abandon),
    );
    pb1.set_prefix("Reading koutput");
    pb1.set_style(reader_style);

    let pb2 = progress.add(ProgressBar::no_length().with_finish(ProgressFinish::Abandon));
    pb2.set_prefix("Writing koutput");
    pb2.set_style(writer_style);

    parse::parse_koutput(
        koutput,
        Some(pb1),
        ofile,
        Some(pb2),
        include_sets,
        exclude_aho,
        compression_level,
        batch_size,
        chunk_bytes,
        nqueue,
        threads,
    )
}
