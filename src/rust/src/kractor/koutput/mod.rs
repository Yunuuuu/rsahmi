use aho_corasick::AhoCorasick;
use memchr::{memchr, memmem};
use rustc_hash::FxHashSet as HashSet;

mod io;
mod mmap;

pub use io::reader_kractor_koutput;
pub use mmap::mmap_kractor_koutput;

static KOUTPUT_TAXID_PREFIX_FINDER: std::sync::LazyLock<memmem::Finder> =
    std::sync::LazyLock::new(|| memmem::Finder::new("(taxid"));

static KOUTPUT_TAXID_PREFIX_FINDERREV: std::sync::LazyLock<memmem::FinderRev> =
    std::sync::LazyLock::new(|| memmem::FinderRev::new("(taxid"));

#[allow(dead_code)]
fn kractor_match_aho(matcher: &AhoCorasick, line: &[u8]) -> bool {
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
                return matcher.find(taxid).is_some();
            } else {
                return false;
            }
        }
        field_index += 1;
        field_start += tab_pos + 1;
    }
    false
}

#[allow(dead_code)]
fn kractor_match_hash(taxid_sets: &HashSet<&[u8]>, line: &[u8]) -> bool {
    // Efficient 3rd column parsing
    let mut field_start = 0usize;
    let mut field_index = 0usize;
    while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
        if field_index == 2 {
            // we don't include the last `\t`
            let taxa = &line[field_start .. (field_start + tab_pos)];
            // extract the taxid pattern
            if let Some(start) = KOUTPUT_TAXID_PREFIX_FINDERREV.rfind(taxa) {
                if let Some(pos) = memchr(b')', &taxa[start ..]) {
                    return taxid_sets.contains(&taxa[start ..= start + pos]);
                };
            };
            return false;
        }
        field_start += tab_pos + 1;
        field_index += 1;
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
        reader_kractor_koutput(
            patterns,
            input_path.to_str().unwrap(),
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
        mmap_kractor_koutput(
            patterns,
            input_path.to_str().unwrap(),
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
