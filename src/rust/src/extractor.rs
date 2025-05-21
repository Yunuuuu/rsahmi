use std::collections::HashSet;
use std::io::{self, BufRead};
use std::path::Path;

use extendr_api::prelude::*;
use futures::TryStreamExt;
use noodles_fasta::r#async::io::Writer;
use noodles_fasta::record::{definition, sequence};
use noodles_fasta::Record;
use noodles_fastq::r#async::io::Reader;
use tokio::fs::File;
use tokio::io::{AsyncWriteExt, BufReader, BufWriter};
use tokio::runtime::Builder;

async fn write_matching_reads<P>(
    fq: P,
    ofile: P,
    id_file: P,
    buffersize: usize,
) -> io::Result<()>
where
    P: AsRef<Path>,
{
    // Open input FASTQ
    let in_file = File::open(fq).await?;

    let in_buf = BufReader::with_capacity(buffersize, in_file);
    let mut reader = Reader::new(in_buf);

    // Open output FASTA
    let out_file = File::create(ofile).await?;
    let out_buf = BufWriter::with_capacity(buffersize, out_file);
    let mut writer = Writer::new(out_buf);

    // Collect IDs into a HashSet for fast lookup
    let id_file = std::fs::File::open(id_file)?;
    let id_buf = std::io::BufReader::with_capacity(buffersize, id_file);
    let id_set = id_buf
        .lines()
        .filter_map(|str| {
            str.ok().and_then(|str| {
                if str.is_empty() {
                    None
                } else {
                    Some(str.as_bytes().to_vec())
                }
            })
        })
        .collect::<HashSet<Vec<u8>>>();

    // Iterate all FASTQ records
    while let Some(record) = reader.records().try_next().await? {
        if id_set.contains(<&[u8]>::from(record.name())) {
            // convert FASTQ to FASTA record and write
            let definition = definition::Definition::new(
                record.name(),
                Some(record.description().to_owned()),
            );
            let sequence =
                sequence::Sequence::from(record.sequence().to_owned());
            let fasta_record = Record::new(definition, sequence);
            writer.write_record(&fasta_record).await?;
        }
    }
    writer.into_inner().flush().await?;
    Ok(())
}

#[extendr]
fn extract_matching_sequence(
    fq: &str,
    ofile: &str,
    id_file: &str,
    buffersize: usize,
) -> io::Result<()> {
    let rt = Builder::new_current_thread().build()?;
    rt.block_on(write_matching_reads(fq, ofile, id_file, buffersize))
}

extendr_module! {
    mod extractor;
    fn extract_matching_sequence;
}
