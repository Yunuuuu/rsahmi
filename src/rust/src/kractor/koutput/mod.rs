use std::fs::File;

use aho_corasick::{AhoCorasick, AhoCorasickKind};
use anyhow::{anyhow, Result};
use indicatif::{ProgressBar, ProgressFinish};
use memchr::{memchr, memchr2, memmem};
use memmap2::Mmap;
use rustc_hash::FxHashMap as HashMap;
use rustc_hash::FxHashSet as HashSet;

pub(crate) mod io;
pub(crate) mod mmap;
use io::reader_kractor_koutput;
use mmap::mmap_kractor_koutput;

use crate::reader::bytes::ProgressBarReader;
use crate::reader::slice::SliceProgressBarReader;

#[allow(clippy::too_many_arguments)]
pub(crate) fn kractor_koutput(
    kreport: &str,
    koutput: &str,
    taxonomy: Option<Vec<&str>>,
    ranks: Option<Vec<&str>>,
    taxa: Option<Vec<&str>>,
    taxids: Option<Vec<&str>>,
    exclude: Option<Vec<&str>>,
    descendants: bool,
    ofile: &str,
    chunk_size: usize,
    buffer_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
    mmap: bool,
) -> Result<()> {
    let mut kreports = crate::kreport::parse_kreport(kreport)?;

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

    let patterns = targeted_taxids
        .into_iter()
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

    let file = File::open(koutput)?;
    if mmap {
        let map = unsafe { Mmap::map(&file) }?;
        crate::mmap_advice(&map)?;
        let reader = SliceProgressBarReader::new(&map);

        mmap_kractor_koutput(
            include_aho,
            exclude_aho,
            reader,
            ofile,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        )
    } else {
        let style = crate::progress_style()?;
        let pb =
            ProgressBar::new(file.metadata()?.len() as u64).with_finish(ProgressFinish::Abandon);
        pb.set_prefix("Parsing koutput");
        pb.set_style(style);
        reader_kractor_koutput(
            include_aho,
            exclude_aho,
            ProgressBarReader::new(&file, pb),
            ofile,
            chunk_size,
            buffer_size,
            batch_size,
            nqueue,
        )
    }
}

pub(crate) static KOUTPUT_TAXID_PREFIX_FINDER: std::sync::LazyLock<memmem::Finder> =
    std::sync::LazyLock::new(|| memmem::Finder::new("(taxid"));

fn kractor_match_aho(
    include_aho: &AhoCorasick,
    exclude_aho: &Option<AhoCorasick>,
    line: &[u8],
) -> bool {
    // println!("Matching line: {:?}", String::from_utf8_lossy(line));
    // Efficient 3rd column parsing
    let mut field_start = 0usize;
    let mut field_index = 0usize;
    while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
        if field_index == 2 {
            // we don't include the last `\t`
            let taxid = &line[field_start .. (field_start + tab_pos)];
            if let Some(start) = KOUTPUT_TAXID_PREFIX_FINDER.find(taxid) {
                let mut input = aho_corasick::Input::new(taxid);
                input.set_start(start);
                if include_aho.find(taxid).is_none() {
                    return false;
                } else if exclude_aho.is_none() {
                    return true;
                }
            } else {
                return false;
            }
        } else if field_index == 3 {
            // Field 4 (LCA): remainder after the last tab
            field_start += tab_pos + 1;
            let lca;
            if let Some(pos) = memchr2(b'\n', b'\t', &line[field_start ..]) {
                lca = &line[field_start .. (field_start + pos)]
            } else {
                lca = &line[field_start ..]
            };
            if let Some(ref exclude_matcher) = exclude_aho {
                return exclude_matcher.find(lca).is_none();
            }
        }
        field_index += 1;
        field_start += tab_pos + 1;
    }
    false
}

#[cfg(test)]
mod tests {
    use std::fs::{read_to_string, File};
    use std::io::Write;

    use anyhow::Result;
    use tempfile::tempdir;

    use super::*;

    fn write_sample_koutput(path: &std::path::Path) -> Result<()> {
        let mut file = File::create(path)?;
        writeln!(file, "read1\tumi1\t(taxid12345)\tother")?;
        writeln!(file, "read2\tumi2\t(taxid54321)\tother")?;
        writeln!(file, "read3\tumi3\t(no_taxid)\tother")?;
        writeln!(file, "read4\tumi4\t(taxid99999)\tother")?;
        Ok(())
    }

    #[test]
    fn test_reader_kractor_koutput() -> Result<()> {
        let dir = tempdir()?;
        let input_path = dir.path().join("input.txt");
        let output_path = dir.path().join("output.txt");

        write_sample_koutput(&input_path)?;

        let patterns = &["taxid12345", "taxid99999"];
        let matcher = AhoCorasick::builder()
            .kind(Some(AhoCorasickKind::DFA))
            .build(patterns)?;
        let file = File::open(input_path).unwrap();
        reader_kractor_koutput(
            matcher,
            None,
            file,
            output_path.to_str().unwrap(),
            2,
            64,
            64,
            Some(10),
        )?;

        let output = read_to_string(output_path)?;
        assert!(output.contains("taxid12345"));
        assert!(output.contains("taxid99999"));
        assert!(!output.contains("taxid54321"));
        assert!(!output.contains("no_taxid"));
        Ok(())
    }

    #[test]
    fn test_mmap_kractor_koutput() -> Result<()> {
        let dir = tempdir()?;
        let input_path = dir.path().join("input.txt");
        let output_path = dir.path().join("output.txt");

        write_sample_koutput(&input_path)?;

        let patterns = &["taxid54321"];
        let matcher = AhoCorasick::builder()
            .kind(Some(AhoCorasickKind::DFA))
            .build(patterns)?;
        let file = File::open(input_path).unwrap();
        let map = unsafe { Mmap::map(&file) }?;
        crate::mmap_advice(&map)?;
        let reader = SliceProgressBarReader::new(&map);
        mmap_kractor_koutput(
            matcher,
            None,
            reader,
            output_path.to_str().unwrap(),
            2,
            64,
            10,
            Some(10),
        )?;

        let output = read_to_string(output_path)?;
        assert!(output.contains("taxid54321"));
        assert!(!output.contains("taxid12345"));
        assert!(!output.contains("taxid99999"));
        Ok(())
    }
}
