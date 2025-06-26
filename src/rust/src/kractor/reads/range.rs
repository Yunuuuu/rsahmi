use anyhow::{anyhow, Result};
use extendr_api::prelude::*;

#[derive(Debug)]
pub enum RangeKind {
    RangeFrom(usize),
    RangeTo(usize),
    Range(usize, usize), // Both start and end are provided
}

impl<'a> std::ops::Index<RangeKind> for [u8] {
    type Output = [u8];

    fn index(&self, index: RangeKind) -> &Self::Output {
        match index {
            RangeKind::RangeFrom(start) => {
                &self[start ..] // Slice from start to the end
            }
            RangeKind::RangeTo(end) => {
                &self[.. end] // Slice from the start to end
            }
            RangeKind::Range(start, end) => {
                &self[start .. end] // Slice from start to end
            }
        }
    }
}

pub fn ubpatterns(patterns: Robj) -> Option<Vec<RangeKind>> {
    patterns.as_list().map_or(None, |r| {
        // Collect valid ranges from the list elements
        let out = r
            .values()
            // for each patten in the Robj
            .filter_map(|pattern| {
                pattern.as_list().map_or(None, |pair| {
                    if pair.len() == 2 {
                        let start = unsafe { pair.get_unchecked(0) };
                        let start_usize = start.as_integer_slice().and_then(|s| {
                            if s.len() == 1 && s[0] >= 0 {
                                Some(s[0] as usize)
                            } else {
                                None
                            }
                        });
                        let end = unsafe { pair.get_unchecked(1) };
                        let end_usize = end.as_integer_slice().and_then(|s| {
                            if s.len() == 1 && s[0] >= 0 {
                                Some(s[0] as usize)
                            } else {
                                None
                            }
                        });
                        if start.is_null() {
                            if let Some(e) = end_usize {
                                Some(RangeKind::RangeTo(e))
                            } else {
                                None
                            }
                        } else if end.is_null() {
                            if let Some(s) = start_usize {
                                Some(RangeKind::RangeFrom(s))
                            } else {
                                None
                            }
                        } else {
                            if let (Some(s), Some(e)) = (start_usize, end_usize) {
                                Some(RangeKind::Range(s, e))
                            } else {
                                None
                            }
                        }
                    } else {
                        None
                    }
                })
            })
            .collect::<Vec<_>>();

        // If no valid ranges were found, return None, otherwise return the list of ranges
        if out.is_empty() {
            None
        } else {
            Some(out)
        }
    })
}

// Function to extract sequence based on patterns and return Result
pub fn extract_pattern_from_sequence(seq: &[u8], ranges: &[RangeKind]) -> Result<Vec<u8>> {
    let mut extracted = Vec::with_capacity(ranges.len());

    for range in ranges {
        match range {
            // Handle RangeFrom (start..)
            RangeKind::RangeFrom(start) => {
                if *start < seq.len() {
                    extracted.push(&seq[*start ..]); // Extract from start to end of sequence
                } else {
                    return Err(anyhow!(
                        "Range {:?} is out of bounds for sequence of length {}",
                        range,
                        seq.len()
                    ));
                }
            }

            // Handle RangeTo (..end)
            RangeKind::RangeTo(end) => {
                if *end <= seq.len() {
                    extracted.push(&seq[.. *end]); // Extract from beginning to end
                } else {
                    return Err(anyhow!(
                        "Range {:?} is out of bounds for sequence of length {}",
                        range,
                        seq.len()
                    ));
                }
            }

            // Handle Range (start..end)
            RangeKind::Range(start, end) => {
                if *start < seq.len() && *end <= seq.len() {
                    extracted.push(&seq[*start .. *end]); // Extract from start to end
                } else {
                    return Err(anyhow!(
                        "Range {:?} is out of bounds for sequence of length {}",
                        range,
                        seq.len()
                    ));
                }
            }
        }
    }

    // Combine all the slices into a single Vec<u8>
    Ok(extracted.concat()) // Using concat() to avoid cloning and flattening
}
