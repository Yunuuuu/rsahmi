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
fn extract_matching_sequence(
    fq1: &str,
    ofile1: &str,
    fq2: Option<&str>,
    ofile2: Option<&str>,
    koutput_file: &str,
    buffersize: usize,
) -> io::Result<()> {
    write_matching_reads(fq1, ofile1, fq2, ofile2, koutput_file, buffersize)
}

fn write_matching_reads<P>(
    fq1: P,
    ofile1: P,
    fq2: Option<P>,
    ofile2: Option<P>,
    koutput_file: P,
    buffersize: usize,
) -> io::Result<()>
where
    P: AsRef<Path> + Display,
{
    // Collect IDs into a HashSet for fast lookup
    rprintln!("Extracting sequence IDs from {}", koutput_file);
    let id_set = read_sequence_id_from_koutput(koutput_file, buffersize)?;

    // Open input FASTQ
    rprintln!("Extracting matching sequence from {}", fq1);
    let in_file1 = File::open(fq1)?;
    let in_buf1 = BufReader::with_capacity(buffersize, in_file1);
    let mut reader1 = Reader::new(in_buf1);

    // Open output FASTA
    let out_file1 = File::create(ofile1)?;
    let out_buf1 = BufWriter::with_capacity(buffersize, out_file1);
    let mut writer1 = Writer::new(out_buf1);

    // Iterate all FASTQ records
    write_matching_records(&mut reader1, &mut writer1, &id_set)?;

    if let (Some(in_file), Some(out_file)) = (fq2, ofile2) {
        rprintln!("Extracting matching sequence from {}", in_file);
        // Open input FASTQ
        let in_file2 = File::open(in_file)?;
        let in_buf2 = BufReader::with_capacity(buffersize, in_file2);
        let mut reader2 = Reader::new(in_buf2);

        // Open output FASTA
        let out_file2 = File::create(out_file)?;
        let out_buf2 = BufWriter::with_capacity(buffersize, out_file2);
        let mut writer2 = Writer::new(out_buf2);
        write_matching_records(&mut reader2, &mut writer2, &id_set)?;
    }
    Ok(())
}

fn read_sequence_id_from_koutput<P>(
    file: P,
    buffersize: usize,
) -> io::Result<HashSet<Vec<u8>>>
where
    P: AsRef<Path>,
{
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

fn write_matching_records<R, W>(
    reader: &mut noodles_fastq::Reader<R>,
    writer: &mut noodles_fasta::Writer<W>,
    id_set: &HashSet<Vec<u8>>,
) -> io::Result<()>
where
    R: BufRead,
    W: Write,
{
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
    mod extractor;
    fn extract_matching_sequence;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::collections::HashSet;
    use tempfile::NamedTempFile;

    #[test]
    fn test_read_sequence_id_from_koutput() -> io::Result<()> {
        // Create a temporary file with tab-separated values
        let mut tmpfile = NamedTempFile::new()?;
        writeln!(tmpfile, "r1\tSEQ001\tdata1")?;
        writeln!(tmpfile, "r2\t\tdata2")?;
        writeln!(tmpfile, "r3\tSEQ002\tdata3")?;
        writeln!(tmpfile, "r4\tSEQ003\tdata4")?;
        writeln!(tmpfile, "r5\t\tdata5")?;
        writeln!(tmpfile, "r6\tSEQ001\tduplicate")?; // Duplicate value

        // Read the file back using your function
        let id_set = read_sequence_id_from_koutput(tmpfile.path(), 1024)?;

        // Define expected results
        let expected: HashSet<Vec<u8>> = vec![
            b"SEQ001".to_vec(),
            b"SEQ002".to_vec(),
            b"SEQ003".to_vec(),
        ]
        .into_iter()
        .collect();

        assert_eq!(id_set, expected);

        Ok(())
    }
}
