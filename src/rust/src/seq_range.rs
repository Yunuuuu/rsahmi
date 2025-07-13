use anyhow::{anyhow, Error, Result};
use extendr_api::prelude::*;

/// A range specification for extracting subsequences from a slice.
///
/// Notes:
/// - Indexing is **1-based** for users, but converted to **0-based** internally.
/// - The `start` is **inclusive**, and `end` is **exclusive**, matching Rust slicing conventions.
/// - This makes the representation consistent with biological conventions while avoiding off-by-one errors.
///
/// Variants:
/// - `RangeFrom(start)`: includes all elements starting from `start` (inclusive).
/// - `RangeTo(end)`: includes all elements up to but not including `end`.
/// - `Range(start, end)`: includes elements from `start` (inclusive) to `end` (exclusive).
#[derive(Debug, Clone)]
pub(crate) enum SeqRange {
    From(usize),
    To(usize),
    Span(usize, usize),
}

impl SeqRange {
    #[allow(dead_code)]
    pub(crate) fn new(start: Option<usize>, end: Option<usize>) -> Self {
        Self::build(start, end).unwrap_or_else(|e| panic!("Failed to create range: {}", e))
    }

    #[allow(dead_code)]
    pub(crate) fn build(start: Option<usize>, end: Option<usize>) -> Result<Self> {
        match (start, end) {
            (Some(s), Some(e)) => {
                if s >= e {
                    Err(anyhow!("start {} is greater than or equal to end {}", s, e))
                } else {
                    Ok(Self::Span(s, e))
                }
            }
            (None, Some(e)) => Ok(Self::To(e)),
            (Some(s), None) => Ok(Self::From(s)),
            (None, None) => {
                return Err(anyhow!(
                    "at least one of 'start' or 'end' must be provided."
                ))
            }
        }
    }

    #[allow(dead_code)]
    pub(crate) fn try_extract<'seq>(&self, sequence: &'seq [u8]) -> Result<&'seq [u8]> {
        let len = sequence.len();
        match *self {
            SeqRange::To(end) => {
                if end > len {
                    Err(anyhow!("RangeTo end {} out of bounds (len = {})", end, len))
                } else {
                    Ok(&sequence[.. end])
                }
            }
            SeqRange::Span(start, end) => {
                if start >= end {
                    Err(anyhow!("Range start {} is greater than end {}", start, end))
                } else if end > len {
                    Err(anyhow!("Range end {} out of bounds (len = {})", end, len))
                } else {
                    Ok(&sequence[start .. end])
                }
            }
            SeqRange::From(start) => {
                if start >= len {
                    Err(anyhow!(
                        "RangeFrom start {} out of bounds (len = {})",
                        start,
                        len
                    ))
                } else {
                    Ok(&sequence[start ..])
                }
            }
        }
    }

    #[allow(dead_code)]
    pub(crate) fn extract<'seq>(&self, sequence: &'seq [u8]) -> &'seq [u8] {
        match *self {
            SeqRange::From(start) => &sequence[start ..],
            SeqRange::To(end) => &sequence[.. end],
            SeqRange::Span(start, end) => &sequence[start .. end],
        }
    }

    #[allow(dead_code)]
    pub(crate) unsafe fn extract_unchecked<'seq>(&self, sequence: &'seq [u8]) -> &'seq [u8] {
        unsafe {
            match *self {
                SeqRange::From(start) => sequence.get_unchecked(start ..),
                SeqRange::To(end) => sequence.get_unchecked(.. end),
                SeqRange::Span(start, end) => sequence.get_unchecked(start .. end),
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct SeqRanges(Vec<SeqRange>);

impl SeqRanges {
    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        Self::with_capacity(0)
    }

    #[allow(dead_code)]
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    #[allow(dead_code)]
    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }

    pub(crate) fn iter(&self) -> std::slice::Iter<SeqRange> {
        self.0.iter()
    }

    #[allow(dead_code)]
    pub(crate) fn sort(&mut self) {
        self.0.sort_by(|x, y| match (x, y) {
            (SeqRange::To(end0), SeqRange::To(end1)) => end0.cmp(end1),
            (SeqRange::To(_), _) => std::cmp::Ordering::Less,
            (_, SeqRange::To(_)) => std::cmp::Ordering::Greater,
            (SeqRange::Span(start0, _), SeqRange::Span(start1, _)) => start0.cmp(start1),
            (SeqRange::From(start0), SeqRange::From(start1)) => start0.cmp(start1),
            (SeqRange::From(_), _) => std::cmp::Ordering::Greater,
            (_, SeqRange::From(_)) => std::cmp::Ordering::Less,
        });
    }

    #[allow(dead_code)]
    pub(crate) fn as_slice(&self) -> &[SeqRange] {
        &self.0
    }
}

/// Checks whether a sorted `SeqRanges` contains overlapping or logically conflicting entries.
///
/// # Requirements:
/// - Input must be **sorted** (`ranges.sort()` must be called beforehand).
/// - Only one `RangeTo` is allowed, and it must be the **first** if used.
/// - Only one `RangeFrom` is allowed, and it must be the **last** if used.
/// - Ranges must not overlap in any form.
pub(crate) fn check_overlap(ranges: &SeqRanges) -> Result<()> {
    // Check for overlapping
    for pair in ranges.0.windows(2) {
        match (&pair[0], &pair[1]) {
            // Case 1: RangeTo must not overlap with a following range
            (SeqRange::To(end0), SeqRange::Span(start1, _))
            | (SeqRange::To(end0), SeqRange::From(start1)) => {
                if end0 > start1 {
                    return Err(anyhow!(
                        "Overlapping ranges: RangeTo({}) overlaps with following range starting at {}.",
                        end0, start1
                    ));
                }
            }

            // Case 2: Only one RangeTo is allowed (must appear first if used)
            (_, SeqRange::To(_)) => {
                // Based on the sort logic, the first must be also RangeTo
                // If the next is RangeTo, they'll always overlap, e.g., [0..10] and [0..20]
                // Hence, multiple RangeTo variants are invalid.
                return Err(anyhow!("Overlapping Ranges: Only one RangeTo is allowed"));
            }

            // Case 3: Overlap between two full ranges or a range and RangeFrom
            (SeqRange::Span(_, end0), SeqRange::Span(start1, _))
            | (SeqRange::Span(_, end0), SeqRange::From(start1)) => {
                if end0 > start1 {
                    return Err(anyhow!(
                        "Overlapping ranges: Range ending at {} overlaps with next range starting at {}.",
                        end0, start1
                    ));
                }
            }

            // Case 4: Only one RangeFrom is allowed (must appear last if used)
            (SeqRange::From(_), _) => {
                // RangeFrom covers [start ..], so any following range would overlap
                return Err(anyhow!("Overlapping Ranges: Only one RangeFrom is allowed"));
            }
        }
    }
    Ok(())
}

impl FromIterator<SeqRange> for SeqRanges {
    fn from_iter<T: IntoIterator<Item = SeqRange>>(iter: T) -> Self {
        iter.into_iter().collect::<Vec<SeqRange>>().into()
    }
}

impl From<Vec<SeqRange>> for SeqRanges {
    fn from(value: Vec<SeqRange>) -> Self {
        Self(value)
    }
}

impl IntoIterator for SeqRanges {
    type IntoIter = std::vec::IntoIter<SeqRange>;
    type Item = SeqRange;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a SeqRanges {
    type IntoIter = std::slice::Iter<'a, SeqRange>;
    type Item = &'a SeqRange;

    fn into_iter(self) -> Self::IntoIter {
        (&self.0).into_iter()
    }
}

/// Extract a vector of `RangeKind` from an R object.
/// The R object must inherit from either `rsahmi_seq_range` or `rsahmi_seq_ranges`.
/// Returns an error if the object is not structured correctly or if the ranges are malformed.
impl TryFrom<&Robj> for SeqRanges {
    type Error = Error;
    fn try_from(value: &Robj) -> Result<Self, Self::Error> {
        if value.inherits("rsahmi_seq_range") {
            // Only one single range
            return SeqRange::try_from(value).map(|range| {
                let mut ranges_vec = SeqRanges::with_capacity(1);
                ranges_vec.0.push(range);
                ranges_vec
            });
        }

        if !value.inherits("rsahmi_seq_ranges") {
            return Err(anyhow!(
            "The object does not inherit a valid range class (expected one of: 'rsahmi_seq_range', or 'rsahmi_seq_ranges')."
        ));
        }

        let list = value
            .as_list()
            .ok_or(anyhow!("Failed to extract valid sequence ranges."))?;

        let mut ranges_vec = SeqRanges::with_capacity(list.len());

        // Iterate over list elements to build vector of RangeKind
        for (i, element) in list.values().enumerate() {
            ranges_vec
                .0
                .push(SeqRange::try_from(&element).map_err(|e| anyhow!("Element {}: {}", i, e))?);
        }
        Ok(ranges_vec)
    }
}

impl TryFrom<&Robj> for SeqRange {
    type Error = Error;
    fn try_from(value: &Robj) -> Result<Self, Self::Error> {
        // Each element should itself be a list of 2 elements
        let pair = value
            .as_list()
            .ok_or_else(|| anyhow!("Robj is not a 2-element list."))?;

        if pair.len() != 2 {
            return Err(anyhow!("Robj must contain exactly 2 values (start, end)."));
        }

        // SAFETY: We already checked that the length is 2
        let start = unsafe { pair.get_unchecked(0) };
        let end = unsafe { pair.get_unchecked(1) };

        let start_usize = if start.is_null() {
            None
        } else {
            let s = start
                .as_integer_slice()
                .ok_or_else(|| anyhow!("start must be an integer."))?;
            if s.len() != 1 || s[0] <= 0 {
                return Err(anyhow!("start must be a single positive integer."));
            }
            Some((s[0] as usize) - 1) // 0-based
        };

        let end_usize = if end.is_null() {
            None
        } else {
            let e = end
                .as_integer_slice()
                .ok_or_else(|| anyhow!("end must be an integer."))?;
            if e.len() != 1 || e[0] <= 0 {
                return Err(anyhow!("end must be a single positive integer."));
            }
            Some(e[0] as usize) // 1-based exclusive
        };
        SeqRange::build(start_usize, end_usize)
    }
}
