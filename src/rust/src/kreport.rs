use std::{fs::File, path::Path};

use anyhow::{anyhow, Context, Result};
use extendr_api::prelude::*;

use crate::{reader0::LineReader, utils::BUFFER_SIZE};

pub(crate) fn parse_kreport<P: AsRef<Path> + ?Sized>(kreport: &P) -> Result<Vec<Kreport>> {
    let path: &Path = kreport.as_ref();
    let mut reader = LineReader::with_capacity(
        BUFFER_SIZE,
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?,
    );
    let mut kreports: Vec<Kreport> = Vec::with_capacity(10);
    let mut ancestors = Vec::with_capacity(10);
    let mut pos = 0; // The line offset of the ancestors
    while let Some(line) = reader.read_line()? {
        let line = line.freeze();
        let fields: Vec<&[u8]> = line.split(|b| *b == b'\t').collect();
        // Parse fixed columns
        let percents = parse_f64(fields[0])?;
        let total_reads = parse_usize(fields[1])?;
        let reads = parse_usize(fields[2])?;
        let minimizer_len;
        let minimizer_n_unique;
        let rank;
        let taxid;
        let taxon;
        let level: usize;
        let mut taxon_field;
        // https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
        // 1. Percentage of fragments covered by the clade rooted at this taxon
        // 2. Number of fragments covered by the clade rooted at this taxon
        // 3. Number of fragments assigned directly to this taxon
        // * 4. Number of minimizers in read data associated with this taxon (new)
        // * 5. An estimate of the number of distinct minimizers in read data
        //    associated with this taxon (new)
        // 6. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
        //    (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that
        //    are not at any of these 10 ranks have a rank code that is formed by
        //    using the rank code of the closest ancestor rank with a number
        //    indicating the distance from that rank. E.g., "G2" is a rank code
        //    indicating a taxon is between genus and species and the grandparent
        //    taxon is at the genus rank.
        // 7. NCBI taxonomic ID number
        // 8. Indented scientific name
        if fields.len() == 6 {
            // 6-column format
            rank = fields[3];
            if rank[0] == b'U' {
                continue;
            }
            taxid = fields[4];
            taxon_field = fields[5].into_iter().peekable();
            minimizer_len = None;
            minimizer_n_unique = None;
        } else if fields.len() == 8 {
            // 8-column format
            rank = fields[5];
            if rank[0] == b'U' {
                continue;
            }
            minimizer_len = Some(parse_usize(fields[3])?);
            minimizer_n_unique = Some(parse_usize(fields[4])?);
            taxid = fields[6];
            taxon_field = fields[7].into_iter().peekable();
        } else {
            return Err(anyhow!(
                "Invalid line with {} fields: {:?}",
                fields.len(),
                String::from_utf8_lossy(&line)
            ))
            .with_context(|| format!("Failed to parse kraken report from {}", path.display()));
        };
        let mut n = 0;
        while let Some(byte) = taxon_field.peek() {
            if **byte == b' ' {
                n += 1;
                taxon_field.next();
            } else {
                break;
            }
        }
        level = n / 2;
        taxon = taxon_field.copied().collect::<Vec<u8>>();
        let rank: Vec<u8> = rank.into_iter().copied().collect();
        let taxid: Vec<u8> = taxid.into_iter().copied().collect();
        while let Some(ancestor) = ancestors.last() {
            if unsafe { kreports.get_unchecked::<usize>(*ancestor) }.level != level - 1 {
                ancestors.pop();
            } else {
                break;
            }
        }
        let ((ranks, taxids), taxa) = ancestors
            .iter()
            .map(|i| {
                let report = unsafe { kreports.get_unchecked::<usize>(*i) };
                (
                    (report.rank.clone(), report.taxid.clone()),
                    report.taxon.clone(),
                )
            })
            .chain(std::iter::once((
                (rank.clone(), taxid.clone()),
                taxon.clone(),
            )))
            .unzip();

        // always remove root species from ancestors
        if rank[0] != b'R' {
            ancestors.push(pos);
        }
        let report = Kreport {
            percents,
            total_reads,
            reads,
            minimizer_len,
            minimizer_n_unique,
            rank,
            taxid,
            taxon,
            ranks,
            taxids,
            taxa,
            level,
        };
        kreports.push(report);
        pos += 1;
    }
    Ok(kreports)
}

#[allow(dead_code)]
pub(crate) struct Kreport {
    pub(crate) percents: f64,
    pub(crate) total_reads: usize,
    pub(crate) reads: usize,
    pub(crate) minimizer_len: Option<usize>,
    pub(crate) minimizer_n_unique: Option<usize>,
    pub(crate) rank: Vec<u8>,
    pub(crate) taxid: Vec<u8>,
    pub(crate) taxon: Vec<u8>,
    pub(crate) ranks: Vec<Vec<u8>>,
    pub(crate) taxids: Vec<Vec<u8>>,
    pub(crate) taxa: Vec<Vec<u8>>,
    pub(crate) level: usize,
}

// Parse &[u8] slice to f64 assuming ASCII decimal representation
fn parse_f64(bytes: &[u8]) -> Result<f64> {
    let s = str::from_utf8(bytes).with_context(|| format!("Invalid UTF-8: {:?}", bytes))?;
    s.trim()
        .parse::<f64>()
        .with_context(|| format!("Failed to parse float '{}'", s))
}

// Parse &[u8] slice to usize assuming ASCII decimal representation
fn parse_usize(bytes: &[u8]) -> Result<usize> {
    let s = str::from_utf8(bytes).with_context(|| format!("Invalid UTF-8: {:?}", bytes))?;
    s.trim()
        .parse::<usize>()
        .with_context(|| format!("Failed to parse integer '{}'", s))
}

fn u8_to_list_rstr(vv: Vec<Vec<u8>>) -> Vec<Rstr> {
    vv.into_iter().map(|v| u8_to_rstr(v)).collect()
}

fn u8_to_rstr(bytes: Vec<u8>) -> Rstr {
    Rstr::from_string(&unsafe { String::from_utf8_unchecked(bytes) })
}

#[extendr]
fn kraken_report(kreport: &str) -> std::result::Result<List, String> {
    let kreports = parse_kreport(kreport).map_err(|e| format!("{:?}", e))?;

    let mut percents = Vec::with_capacity(kreports.len());
    let mut total_reads = Vec::with_capacity(kreports.len());
    let mut reads = Vec::with_capacity(kreports.len());

    let mut rank = Vec::with_capacity(kreports.len());
    let mut taxid = Vec::with_capacity(kreports.len());
    let mut taxa = Vec::with_capacity(kreports.len());

    let mut ranks = Vec::with_capacity(kreports.len());
    let mut taxids = Vec::with_capacity(kreports.len());
    let mut taxon = Vec::with_capacity(kreports.len());

    // Optional columns
    let mut minimizer_len = Vec::with_capacity(kreports.len());
    let mut minimizer_n_unique = Vec::with_capacity(kreports.len());
    for report in kreports {
        percents.push(report.percents);
        total_reads.push(report.total_reads as f64);
        reads.push(report.reads as f64);

        rank.push(u8_to_rstr(report.rank));
        taxid.push(u8_to_rstr(report.taxid));
        taxa.push(u8_to_rstr(report.taxon));

        ranks.push(Robj::from(u8_to_list_rstr(report.ranks)));
        taxids.push(Robj::from(u8_to_list_rstr(report.taxids)));
        taxon.push(Robj::from(u8_to_list_rstr(report.taxa)));

        match (report.minimizer_len, report.minimizer_n_unique) {
            (Some(a), Some(b)) => {
                minimizer_len.push(a as f64);
                minimizer_n_unique.push(b as f64);
            }
            _ => {}
        }
    }

    let ranks = List::from_values(ranks);
    let taxids = List::from_values(taxids);
    let taxon = List::from_values(taxon);

    // Create R dataframe
    let out = if minimizer_len.is_empty() {
        list![
            percents = percents,
            total_reads = total_reads,
            reads = reads,
            rank = rank,
            taxid = taxid,
            taxa = taxa,
            ranks = ranks,
            taxids = taxids,
            taxon = taxon
        ]
    } else {
        list![
            percents = percents,
            total_reads = total_reads,
            reads = reads,
            minimizer_len = minimizer_len,
            minimizer_n_unique = minimizer_n_unique,
            rank = rank,
            taxid = taxid,
            taxa = taxa,
            ranks = ranks,
            taxids = taxids,
            taxon = taxon
        ]
    };
    Ok(out)
}

extendr_module! {
    mod kreport;
    fn kraken_report;
}
