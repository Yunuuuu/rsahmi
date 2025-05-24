use aho_corasick::AhoCorasick;
use memchr::memchr;

use crate::chunk::{ChunkParser, ChunkProcessor};

pub struct KOutputProcessor {
    matcher: AhoCorasick,
}

impl KOutputProcessor {
    pub fn new(matcher: AhoCorasick) -> Self {
        Self { matcher }
    }
}

pub struct KOutputParser<'a> {
    matcher: &'a AhoCorasick,
}

impl<'a> ChunkProcessor for &'a KOutputProcessor {
    type Parser = KOutputParser<'a>;

    fn new_parser(&self) -> Self::Parser {
        KOutputParser::new(&self.matcher)
    }
}

impl ChunkParser for KOutputParser<'_> {
    fn parse<F>(
        &self,
        chunk: Vec<u8>,
        push: &mut F,
    ) -> std::result::Result<(), String>
    where
        F: FnMut(Vec<u8>) -> std::result::Result<(), String>,
    {
        let mut start = 0;
        while let Some(line_pos) = memchr(b'\n', &chunk[start ..]) {
            // we include the last `\n` for writing
            let line = &chunk[start ..= (start + line_pos)];
            self.parse_line(line, push)?;
            // push(line.to_vec())?;
            start += line_pos + 1;
        }
        Ok(())
    }
}

impl<'a> KOutputParser<'a> {
    fn new(matcher: &'a AhoCorasick) -> Self {
        Self { matcher }
    }

    fn parse_line<F>(
        &self,
        line: &[u8],
        push: &mut F,
    ) -> std::result::Result<(), String>
    where
        F: FnMut(Vec<u8>) -> std::result::Result<(), String>,
    {
        // Efficient 3rd column parsing
        let mut field_start = 0usize;
        let mut field_count = 0usize;
        while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
            if field_count == 2 {
                // we don't include the last `\t`
                let taxid = &line[field_start .. (field_start + tab_pos)];
                if self.matcher.find(taxid).is_some() {
                    push(line.to_vec())?;
                }
                break;
            }
            field_start += tab_pos + 1;
            field_count += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use aho_corasick::AhoCorasick;

    use super::*;

    #[test]
    fn test_koutput_parser_filters_by_third_column() {
        let matcher = AhoCorasick::new(["123"]).unwrap();
        let processor = KOutputProcessor::new(matcher);
        let parser = (&processor).new_parser();

        let input = b"id1\tname\t123\textra\nid2\tname\t999\textra\n".to_vec();

        let mut output = Vec::new();
        parser
            .parse(input, &mut |line| {
                output.push(line);
                Ok(())
            })
            .unwrap();
        assert_eq!(output.len(), 1);
        assert_eq!(output[0], b"id1\tname\t123\textra\n".to_vec());
    }
}
