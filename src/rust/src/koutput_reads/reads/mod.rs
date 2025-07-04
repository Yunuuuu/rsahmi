mod io;
mod mmap;
use rustc_hash::FxHashMap as HashMap;

/// Return `true` if all base counts are ≤ `threshold`, otherwise `false`.
fn pass_complexity_filter(seq: &[u8], threshold: usize) -> bool {
    // remove low complexity reads (<20 non-sequentially repeated nucleotides)
    let threshold = seq.len() - threshold;
    let mut counts = HashMap::with_capacity_and_hasher(4, rustc_hash::FxBuildHasher); // ATGC
    for &b in seq {
        let count = counts.entry(b).or_insert(0);
        *count += 1;
        if *count > threshold {
            return false; // Early exit: too many repeats
        }
    }
    true
}

/// Returns `true` if all quality scores are ≥ `min_phred`.
fn pass_quality_filter(qual: &[u8], threshold: u8) -> bool {
    // threshold 53 for Phred score < 20 (Phred+33 ASCII)
    // threshold 84 for Phred score < 20 (Phred+64 ASCII)
    qual.iter().all(|&q| q >= threshold)
}
