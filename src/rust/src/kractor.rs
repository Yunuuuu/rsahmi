use std::collections::HashSet;
use std::fmt::Display;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use extendr_api::prelude::*;
use noodles_fasta::io::Writer;
use noodles_fasta::record::{definition, sequence, Record};
use noodles_fastq::io::Reader;

#[extendr]
#[allow(clippy::too_many_arguments)]
fn kractor(
    koutput: &str,
    taxids: Robj,
    ofile: &str,
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    buffersize: usize,
) -> std::result::Result<(), String> {
    let taxids: Vec<&str> = taxids
        .as_str_vector()
        .ok_or("`taxids` must be a character vector")?;
    let id_set = write_matching_output(koutput, &taxids, ofile, buffersize)
        .map_err(|e| e.to_string())?;
    write_matching_reads(fq1, ofile1, fq2, ofile2, &id_set, buffersize)
        .map_err(|e| e.to_string())
}

fn write_matching_output<P>(
    koutput: P,
    taxids: &[&str],
    ofile: P,
    buffersize: usize,
) -> io::Result<HashSet<Vec<u8>>>
where
    P: AsRef<Path> + Display,
{
    rprintln!("Extracting matching kraken2 output {}", koutput);
    let input = BufReader::with_capacity(buffersize, File::open(koutput)?);
    let mut output = BufWriter::with_capacity(buffersize, File::create(ofile)?);

    // Preallocate with expected number of sequence IDs
    let id_set: HashSet<Vec<u8>> = input
        .lines()
        .filter_map(|read| {
            match read {
                Ok(line) => {
                    let mut fields = line.split('\t');
                    fields.next(); // Just omit the first item
                    let second = fields.next(); // sequence ID
                    let third = fields.next(); // taxids
                    match third {
                        Some(taxid)
                            if taxids.iter().any(|p| taxid.contains(p)) =>
                        {
                            // write the matching output
                            let _ = writeln!(output, "{line}");

                            // extract the sequence id
                            second.and_then(|seq_id| {
                                // we remove empty sequence IDs
                                if !seq_id.is_empty() {
                                    Some(seq_id.as_bytes().to_vec())
                                } else {
                                    None
                                }
                            })
                        }
                        _ => None,
                    }
                }
                _ => None,
            }
        })
        .collect::<HashSet<Vec<u8>>>();
    output.flush()?;
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
        rprintln!("Extracting matching sequence from {}", in_file);
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
