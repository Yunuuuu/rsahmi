use std::collections::HashSet;
use std::fmt::Display;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use extendr_api::prelude::*;
use memchr::memmem;
use noodles_fasta::io::Writer;
use noodles_fasta::record::{definition, sequence, Record};
use noodles_fastq::io::Reader;

#[extendr]
#[allow(clippy::too_many_arguments)]
fn kractor(
    koutput: &str,
    patterns: Robj,
    ofile: &str,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    buffersize: usize,
) -> std::result::Result<(), String> {
    let patterns: Vec<&str> = patterns
        .as_str_vector()
        .ok_or("`patterns` must be a character vector")?;
    write_matching_output(koutput, &patterns, ofile, buffersize)
        .map_err(|e| e.to_string())?;
    let id_set = read_sequence_id_from_koutput(ofile, buffersize)
        .map_err(|e| e.to_string())?;
    write_matching_reads(fq1, ofile1, fq2, ofile2, &id_set, buffersize)
        .map_err(|e| e.to_string())
}

fn write_matching_output<P>(
    koutput: P,
    patterns: &[&str],
    ofile: P,
    buffersize: usize,
) -> io::Result<()>
where
    P: AsRef<Path> + Display,
{
    rprintln!("Extracting matching kraken2 output from: {}", koutput);
    let mut input = BufReader::with_capacity(buffersize, File::open(koutput)?);
    let mut output = BufWriter::with_capacity(buffersize, File::create(ofile)?);
    let needles: Vec<_> = patterns
        .iter()
        .map(|s| memmem::Finder::new(s.as_bytes()))
        .collect();
    let mut line = Vec::new();
    while let Ok(bytes_read) = input.read_until(b'\n', &mut line) {
        if bytes_read == 0 {
            break;
        }
        let mut fields = line.split(|str| *str == b'\t');
        // second field; // sequence ID
        // third field; // taxids
        if let Some(taxid) = fields.nth(2) {
            if needles.iter().any(|needle| needle.find(taxid).is_some()) {
                output.write_all(&line)?;
            }
        }
        line.clear()
    }
    Ok(())
}

fn read_sequence_id_from_koutput<P>(
    file: P,
    buffersize: usize,
) -> io::Result<HashSet<Vec<u8>>>
where
    P: AsRef<Path> + Display,
{
    // Collect IDs into a HashSet for fast lookup
    rprintln!("Extracting sequence IDs from: {}", file);
    let opened = File::open(file)?;
    let buffer = BufReader::with_capacity(buffersize, opened);
    let id_set = buffer
        .lines()
        .filter_map(|line| {
            line.ok().and_then(|str| {
                // we selected the second column
                str.split("\t").nth(1).and_then(|second| {
                    // we remove empty sequence IDs
                    if second.is_empty() {
                        None
                    } else {
                        Some(second.as_bytes().to_vec())
                    }
                })
            })
        })
        .collect::<HashSet<Vec<u8>>>();
    Ok(id_set)
}

fn write_matching_reads<P>(
    fq1: P,
    ofile1: P,
    fq2: Option<P>,
    ofile2: Option<P>,
    id_set: &HashSet<Vec<u8>>,
    buffersize: usize,
) -> io::Result<()>
where
    P: AsRef<Path> + Display,
{
    // Iterate all FASTQ records
    write_matching_records(fq1, ofile1, buffersize, id_set)?;

    if let (Some(in_file), Some(out_file)) = (fq2, ofile2) {
        rprintln!("Extracting the matching sequence from: {}", in_file);
        write_matching_records(in_file, out_file, buffersize, id_set)?;
    }
    Ok(())
}

fn write_matching_records<P>(
    fastq: P,
    fasta: P,
    buffersize: usize,
    id_set: &HashSet<Vec<u8>>,
) -> io::Result<()>
where
    P: AsRef<Path> + Display,
{
    // Open input FASTQ
    let in_file = File::open(fastq)?;
    let in_buf = BufReader::with_capacity(buffersize, in_file);
    let mut reader = Reader::new(in_buf);

    // Open output FASTA
    let out_file = File::create(fasta)?;
    let out_buf = BufWriter::with_capacity(buffersize, out_file);
    let mut writer = Writer::new(out_buf);

    // Iterate all FASTQ records
    for read in reader.records() {
        let record = read?;
        if id_set.contains(<&[u8]>::from(record.name())) {
            // convert FASTQ to FASTA record and write
            let definition = definition::Definition::new(
                record.name(),
                Some(record.description().to_owned()),
            );
            let sequence =
                sequence::Sequence::from(record.sequence().to_owned());
            let fasta_record = Record::new(definition, sequence);
            writer.write_record(&fasta_record)?;
        }
    }
    writer.get_mut().flush()?;
    Ok(())
}

extendr_module! {
    mod kractor;
    fn kractor;
}
