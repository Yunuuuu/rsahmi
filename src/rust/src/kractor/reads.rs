use std::collections::HashSet;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use memchr::memchr;

use crate::chunk::ChunkSplitter;
use crate::chunk::{ChunkParser, ChunkProcessor};

pub struct ReadsProcessor<'a> {
    sequence_ids: &'a HashSet<Vec<u8>>,
}

impl<'a> ReadsProcessor<'a> {
    pub fn new(ids: &'a HashSet<Vec<u8>>) -> Self {
        Self { sequence_ids: ids }
    }
}

pub struct ReadsParser<'a> {
    sequence_ids: &'a HashSet<Vec<u8>>,
}

impl<'a> ChunkProcessor for ReadsProcessor<'a> {
    type Parser = ReadsParser<'a>;

    fn new_parser(&self) -> Self::Parser {
        ReadsParser::new(self.sequence_ids)
    }
    fn new_splitter(&self) -> ChunkSplitter {
        ChunkSplitter::Bytes(b"\n@")
    }
}

#[allow(unused_assignments)]
impl<'a> ChunkParser for ReadsParser<'a> {
    fn parse<F>(
        &self,
        chunk: Vec<u8>,
        push: &mut F,
    ) -> std::result::Result<(), String>
    where
        F: FnMut(Vec<u8>) -> std::result::Result<(), String>,
    {
        let mut start = 0;
        let mut record_pos = 0usize;
        let mut record_start = 0usize;
        let mut sequence_start = 0usize;
        let mut sequence_end = 0usize;
        let mut push_record = false;
        while let Some(line_pos) = memchr(b'\n', &chunk[start ..]) {
            match record_pos {
                0 => {
                    // first line, should be a sequence ID
                    record_start = start;
                    // check sequence is valid
                    if chunk[start] != b'@' {
                        push_record = false;
                    } else {
                        // remove the '@' from the start of the sequence ID
                        // ID is after '@' and before first space
                        // id and description were split by a space, so we take the first part
                        let end = memchr(
                            b' ',
                            &chunk[(start + 1) ..= (start + line_pos)],
                        )
                        .map_or_else(
                            // NO description
                            || start + line_pos,
                            // we don't add 1 here since we don't want to include the space
                            |offset| start + offset,
                        );
                        push_record = self
                            .sequence_ids
                            .contains(&chunk[(start + 1) ..= end]);
                    }
                    record_pos += 1;
                }
                1 => {
                    // second line, should be the sequence
                    record_pos += 1;
                    sequence_start = start;
                    sequence_end = start + line_pos;
                }
                2 => {
                    // third line, should be a plus line
                    if push_record && chunk[start] != b'+' {
                        push_record = false;
                    }
                    record_pos += 1;
                }
                3 => {
                    // fourth line, should be the quality scores
                    // It should have the same length as the sequence
                    // sequence and quality lengths do not match, skip this record
                    if push_record
                        && (sequence_end - sequence_start) == line_pos
                    {
                        // push the id, description, and sequence
                        // we remove the '@' from the start of the sequence ID
                        let record =
                            &chunk[(record_start + 1) ..= sequence_end];
                        push(record.to_vec())?;
                    }
                    record_pos = 0; // reset for next record
                }
                _ => unreachable!(),
            }
            start += line_pos + 1;
        }
        Ok(())
    }
}

impl<'a> ReadsParser<'a> {
    fn new(ids: &'a HashSet<Vec<u8>>) -> Self {
        Self { sequence_ids: ids }
    }
}

pub fn read_sequence_id_from_koutput<P>(
    file: P,
    buffersize: usize,
) -> std::result::Result<HashSet<Vec<u8>>, String>
where
    P: AsRef<Path> + Display,
{
    let opened =
        File::open(file).map_err(|e| format!("Open file failed: {}", e))?;
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

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn test_reads_parser_with_matching_id() {
        let ids: HashSet<Vec<u8>> = [b"seq1".to_vec(), b"seq2".to_vec()]
            .iter()
            .cloned()
            .collect();
        let processor = ReadsProcessor::new(&ids);
        let parser = processor.new_parser();

        // Simulate a complete FASTQ record for seq1 and one unmatched record
        let fastq_data = b"@seq1 description\nACGTACGT\n+\nIIIIIIII\n@seq3\nGATTACA\n+\nIIIIIII\n";

        let mut results = Vec::new();
        parser
            .parse(fastq_data.to_vec(), &mut |record| {
                results.push(record);
                Ok(())
            })
            .unwrap();
        // Only the first record should be collected
        assert_eq!(results.len(), 1);
        assert_eq!(results[0], b"seq1 description\nACGTACGT\n".to_vec());
    }

    #[test]
    fn test_reads_parser_with_invalid_format() {
        let ids: HashSet<Vec<u8>> =
            [b"seq1".to_vec()].iter().cloned().collect();
        let processor = ReadsProcessor::new(&ids);
        let parser = processor.new_parser();

        // Corrupt FASTQ: missing '+' line
        let fastq_data = b"@seq1\nACGTACGT\nIIIIIIII\n";

        let mut results = Vec::new();
        parser
            .parse(fastq_data.to_vec(), &mut |record| {
                results.push(record);
                Ok(())
            })
            .unwrap();

        // Should not collect anything
        assert!(results.is_empty());
    }
}
