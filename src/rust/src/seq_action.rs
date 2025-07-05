use anyhow::{anyhow, Result};
use bytes::{BufMut, Bytes, BytesMut};
use extendr_api::prelude::*;

use crate::parser::fastq::FastqRecord;

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
enum SeqRange {
    From(usize),
    To(usize),
    Span(usize, usize),
}

impl SeqRange {
    fn new_checked(start: Option<usize>, end: Option<usize>) -> std::result::Result<Self, String> {
        match (start, end) {
            (Some(s), Some(e)) => {
                if s >= e {
                    Err(format!("start {} is greater than or equal to end {}", s, e))
                } else {
                    Ok(Self::Span(s, e))
                }
            }
            (None, Some(e)) => Ok(Self::To(e)),
            (Some(s), None) => Ok(Self::From(s)),
            (None, None) => {
                return Err("at least one of 'start' or 'end' must be provided.".to_string())
            }
        }
    }

    fn new(start: Option<usize>, end: Option<usize>) -> Self {
        Self::new_checked(start, end).unwrap_or_else(|e| panic!("Failed to create range: {}", e))
    }

    fn try_extract<'s>(&self, slice: &'s [u8]) -> Result<&'s [u8]> {
        let len = slice.len();
        match *self {
            SeqRange::To(end) => {
                if end > len {
                    Err(anyhow!("RangeTo end {} out of bounds (len = {})", end, len))
                } else {
                    Ok(&slice[.. end])
                }
            }
            SeqRange::Span(start, end) => {
                if start >= end {
                    Err(anyhow!("Range start {} is greater than end {}", start, end))
                } else if end > len {
                    Err(anyhow!("Range end {} out of bounds (len = {})", end, len))
                } else {
                    Ok(&slice[start .. end])
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
                    Ok(&slice[start ..])
                }
            }
        }
    }

    fn extract<'s>(&self, slice: &'s [u8]) -> &'s [u8] {
        match *self {
            SeqRange::From(start) => &slice[start ..],
            SeqRange::To(end) => &slice[.. end],
            SeqRange::Span(start, end) => &slice[start .. end],
        }
    }

    unsafe fn extract_unchecked<'s>(&self, slice: &'s [u8]) -> &'s [u8] {
        unsafe {
            match *self {
                SeqRange::From(start) => slice.get_unchecked(start ..),
                SeqRange::To(end) => slice.get_unchecked(.. end),
                SeqRange::Span(start, end) => slice.get_unchecked(start .. end),
            }
        }
    }
}

#[derive(Debug, Clone)]
struct SeqRanges(Vec<SeqRange>);

impl SeqRanges {
    fn new() -> Self {
        Self::with_capacity(0)
    }

    fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn sort(&mut self) {
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

    fn extend<I: IntoIterator<Item = SeqRange>>(&mut self, iter: I) {
        self.0.extend(iter);
    }

    fn push(&mut self, value: SeqRange) {
        self.0.push(value);
    }

    fn as_slice(&self) -> &[SeqRange] {
        &self.0
    }

    fn into_inner(self) -> Vec<SeqRange> {
        self.0
    }
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

// DOn't allow add new sequence range
#[derive(Debug, Clone)]
struct SortedSeqRanges(Vec<SeqRange>);

impl From<SeqRanges> for SortedSeqRanges {
    fn from(mut value: SeqRanges) -> Self {
        value.sort();
        SortedSeqRanges(value.0)
    }
}

impl FromIterator<SeqRange> for SortedSeqRanges {
    fn from_iter<T: IntoIterator<Item = SeqRange>>(iter: T) -> Self {
        let ranges: SeqRanges = iter.into_iter().collect::<Vec<SeqRange>>().into();
        ranges.into()
    }
}

impl SortedSeqRanges {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn as_slice(&self) -> &[SeqRange] {
        &self.0
    }

    fn check_overlap(&self) -> Result<()> {
        // Check for overlapping
        for pair in self.0.windows(2) {
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
}

impl IntoIterator for SortedSeqRanges {
    type IntoIter = std::vec::IntoIter<SeqRange>;
    type Item = SeqRange;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a SortedSeqRanges {
    type IntoIter = std::slice::Iter<'a, SeqRange>;
    type Item = &'a SeqRange;

    fn into_iter(self) -> Self::IntoIter {
        (&self.0).into_iter()
    }
}

struct SubseqEmbedAction {
    tag: String,
    ranges: SortedSeqRanges,
}

impl SubseqEmbedAction {
    fn new(tag: String, ranges: SortedSeqRanges) -> Self {
        Self { tag, ranges }
    }

    fn tag(&self, seq: &[u8]) -> Result<Bytes> {
        let mut sequences = Vec::with_capacity(self.ranges.len());
        for range in &self.ranges {
            sequences.push(range.try_extract(seq)?);
        }
        Ok(self
            .tag
            .as_bytes()
            .iter()
            // add ":" between tag and extracted sequence
            .chain(std::iter::once(&b':'))
            .chain(sequences.into_iter().flatten())
            .copied()
            .collect())
    }
}

struct SubseqTrimAction {
    ranges: SortedSeqRanges,
}

impl SubseqTrimAction {
    fn new(sorted: SortedSeqRanges) -> Self {
        Self { ranges: sorted }
    }

    /// Trim the sequence and quality string using the defined `RangeKind`s.
    /// If any range is out of bounds, it will return an error.
    fn trim(&self, seq: &[u8], qual: &[u8]) -> Result<(Bytes, Bytes)> {
        let len = seq.len();
        let mut keep = Vec::new();
        let mut cursor = 0;

        // Iterate through the ranges, collecting segments to keep.
        for range in &self.ranges {
            match *range {
                SeqRange::To(end) => {
                    if end > len {
                        return Err(anyhow!("RangeTo end {} out of bounds (len = {})", end, len));
                    } else if cursor < end {
                        cursor = end; // Skip the region up to `end`
                    }
                }
                SeqRange::Span(start, end) => {
                    if start >= end {
                        return Err(anyhow!("Range start {} is greater than end {}", start, end));
                    } else if end > len {
                        return Err(anyhow!("Range end {} out of bounds (len = {})", end, len));
                    } else if start > cursor {
                        keep.push((cursor, start)); // Keep everything before `start`
                        cursor = end; // Skip the region from `start` to `end`
                    } else if cursor < end {
                        cursor = end;
                    }
                }
                SeqRange::From(start) => {
                    if start >= len {
                        return Err(anyhow!(
                            "RangeFrom start {} out of bounds (len = {})",
                            start,
                            len
                        ));
                    } else if start > cursor {
                        keep.push((cursor, start)); // Keep everything before `start`
                    }
                    cursor = len;
                    break;
                }
            }
        }
        if cursor < len {
            keep.push((cursor, len))
        }
        let trimmed_seq = keep
            .iter()
            .map(|(start, end)| &seq[*start .. *end])
            .flatten()
            .copied()
            .collect();
        let trimmed_qual = keep
            .iter()
            .map(|(start, end)| &qual[*start .. *end])
            .flatten()
            .copied()
            .collect();
        Ok((trimmed_seq, trimmed_qual))
    }
}

pub(crate) struct SubseqActions {
    embed: Vec<SubseqEmbedAction>,
    trim: SubseqTrimAction,
}

impl SubseqActions {
    pub(crate) fn transform_fastq_slice(
        &self,
        record: &FastqRecord<&[u8]>,
    ) -> Result<FastqRecord<Bytes>> {
        let description;
        if self.embed.is_empty() {
            description = record.desc.map(|d| Bytes::copy_from_slice(d));
        } else {
            let mut tags = Vec::with_capacity(self.embed.len());
            for action in &self.embed {
                tags.push(action.tag(record.seq)?);
            }

            // prepare the description
            let prefix = b"RSAHMI{";
            let suffix = b'}';
            let mut tag = BytesMut::with_capacity(
                // original description length + prefix length + all tags + separator between each tag + suffix length
                record.desc.map_or(0, |d| d.len() + 1)
                    + prefix.len()
                    + tags.iter().map(|t| t.len()).sum::<usize>()
                    + tags.len()
                    - 1
                    + 1,
            );
            if let Some(desc) = record.desc {
                tag.extend_from_slice(desc);
                tag.put_u8(b' ');
            }
            tag.extend_from_slice(prefix);
            for (i, t) in tags.iter().enumerate() {
                if i > 0 {
                    tag.put_u8(b':');
                }
                tag.extend_from_slice(t);
            }
            tag.put_u8(suffix);
            description = Some(tag.freeze());
        }

        // prepare sequence and quality
        let sequence;
        let quality;
        if self.trim.ranges.len() > 0 {
            (sequence, quality) = self.trim.trim(record.seq, record.qual)?;
        } else {
            sequence = Bytes::copy_from_slice(record.seq);
            quality = Bytes::copy_from_slice(record.qual);
        }

        // construct the output
        Ok(FastqRecord::new(
            Bytes::copy_from_slice(record.id),
            description,
            sequence,
            Bytes::copy_from_slice(record.sep),
            quality,
        ))
    }

    pub(crate) fn transform_fastq_bytes(&self, record: &mut FastqRecord<Bytes>) -> Result<()> {
        if !self.embed.is_empty() {
            let mut tags = Vec::with_capacity(self.embed.len());
            for action in &self.embed {
                tags.push(action.tag(&record.seq)?);
            }

            // prepare the description
            let prefix = b"RSAHMI{";
            let suffix = b'}';
            let mut tag = BytesMut::with_capacity(
                // original description length + prefix length + all tags + separator between each tag + suffix length
                record.desc.as_ref().map_or(0, |d| d.len() + 1)
                    + prefix.len()
                    + tags.iter().map(|t| t.len()).sum::<usize>()
                    + tags.len()
                    - 1
                    + 1,
            );
            if let Some(desc) = &record.desc {
                tag.extend_from_slice(&desc);
                tag.put_u8(b' ');
            }
            tag.extend_from_slice(prefix);
            for (i, t) in tags.iter().enumerate() {
                if i > 0 {
                    tag.put_u8(b':');
                }
                tag.extend_from_slice(t);
            }
            tag.put_u8(suffix);
            record.desc = Some(tag.freeze());
        }

        // prepare sequence and quality
        if self.trim.ranges.len() > 0 {
            let (sequence, quality) = self.trim.trim(&record.seq, &record.qual)?;
            record.seq = sequence;
            record.qual = quality;
        }

        Ok(())
    }
}

struct SubseqActionsBuilder {
    embed_list: Vec<(String, SortedSeqRanges)>,
    trim_list: Vec<SortedSeqRanges>,
}

impl SubseqActionsBuilder {
    fn new() -> Self {
        Self {
            embed_list: Vec::new(),
            trim_list: Vec::new(),
        }
    }

    fn add_embed(&mut self, tag: String, ranges: SortedSeqRanges) {
        self.embed_list.push((tag, ranges));
    }

    fn add_trim(&mut self, ranges: SortedSeqRanges) {
        self.trim_list.push(ranges);
    }

    fn add_action(&mut self, action: SeqAction, ranges: SortedSeqRanges) {
        match action {
            SeqAction::Embed(tag) => {
                self.add_embed(tag, ranges);
            }
            SeqAction::Trim => {
                self.add_trim(ranges);
            }
            SeqAction::EmbedTrim(tag) => {
                self.add_embed(tag, ranges.clone());
                self.add_trim(ranges);
            }
        }
    }

    fn build(self) -> SubseqActions {
        let embed_actions = self
            .embed_list
            .into_iter()
            .map(|(tag, ranges)| SubseqEmbedAction::new(tag, ranges))
            .collect();
        let sorted: SortedSeqRanges = self.trim_list.into_iter().flatten().collect();
        SubseqActions {
            embed: embed_actions,
            trim: SubseqTrimAction::new(sorted),
        }
    }
}

enum SeqAction {
    Embed(String),
    Trim,
    EmbedTrim(String),
}

// Create object from R
pub(crate) fn robj_to_seq_range_actions<'r>(ranges: &Robj, label: &str) -> Result<SubseqActions> {
    let mut out = SubseqActionsBuilder::new();
    for ref robj in ranges
        .as_list()
        .ok_or(anyhow!("Expected a list of sequence range objects."))?
        .values()
    {
        let (action, sorted) = robj_dissect_range_obj(robj, label)?;
        out.add_action(action, sorted)
    }
    Ok(out.build())
}

fn extract_tag(ranges: &Robj) -> &str {
    ranges
        .get_attrib("tag")
        .and_then(|t| t.as_str())
        .unwrap_or_else(|| "rsahmi")
}

fn extract_ordering(ranges: &Robj) -> Option<Vec<usize>> {
    ranges.get_attrib("ordering").and_then(|o| {
        if let Some(ordering) = o.as_integer_slice() {
            if ordering.iter().any(|o| *o < 0) {
                None
            } else {
                Some(ordering.iter().map(|o| *o as usize).collect())
            }
        } else {
            None
        }
    })
}

fn robj_dissect_range_obj(ranges: &Robj, label: &str) -> Result<(SeqAction, SortedSeqRanges)> {
    let seq_ranges = robj_to_ranges(ranges, label)?;

    let sorted: SortedSeqRanges = seq_ranges.into();
    sorted.check_overlap()?;

    let resolved_action = if ranges.inherits("rsahmi_trim") {
        SeqAction::Trim
    } else if ranges.inherits("rsahmi_embed") {
        SeqAction::Embed(extract_tag(ranges).to_string())
    } else if ranges.inherits("rsahmi_embed_trim") {
        SeqAction::EmbedTrim(extract_tag(ranges).to_string())
    } else {
        return Err(anyhow!(
            "The object does not inherit a valid range class (expected one of: 'rsahmi_trim', 'rsahmi_embed', or 'rsahmi_embed_trim')."
        ));
    };
    Ok((resolved_action, sorted))
}

/// Extract a vector of `RangeKind` from an R object.
/// The R object must inherit from either `rsahmi_seq_range` or `rsahmi_seq_ranges`.
/// Returns an error if the object is not structured correctly or if the ranges are malformed.
fn robj_to_ranges(ranges: &Robj, label: &str) -> Result<SeqRanges> {
    if ranges.inherits("rsahmi_seq_range") {
        // Only one single range
        return robj_to_range(ranges, 0, label).map(|range| {
            let mut ranges = SeqRanges::with_capacity(1);
            ranges.push(range);
            ranges
        });
    }

    if !ranges.inherits("rsahmi_seq_ranges") {
        return Err(anyhow!(
            "The object does not inherit a valid range class (expected one of: 'rsahmi_seq_range', or 'rsahmi_seq_ranges')."
        ));
    }

    let list = ranges
        .as_list()
        .ok_or(anyhow!("Failed to extract valid sequence ranges."))?;

    let mut ranges_list = SeqRanges::with_capacity(list.len());

    // Iterate over list elements to build vector of RangeKind
    for (i, element) in list.values().enumerate() {
        ranges_list.push(robj_to_range(&element, i, label)?);
    }
    Ok(ranges_list)
}

fn robj_to_range(range: &Robj, i: usize, label: &str) -> Result<SeqRange> {
    // Each element should itself be a list of 2 elements
    let pair = range
        .as_list()
        .ok_or_else(|| anyhow!("{}: Element {} is not a 2-element list.", label, i + 1))?;

    if pair.len() != 2 {
        return Err(anyhow!(
            "{}: Element {} must contain exactly 2 values (start, end).",
            label,
            i + 1
        ));
    }

    // SAFETY: We already checked that the length is 2
    let start = unsafe { pair.get_unchecked(0) };
    let end = unsafe { pair.get_unchecked(1) };

    let start_usize = if start.is_null() {
        None
    } else {
        let s = start
            .as_integer_slice()
            .ok_or_else(|| anyhow!("{}: Element {} start must be an integer.", label, i + 1))?;
        if s.len() != 1 || s[0] <= 0 {
            return Err(anyhow!(
                "{}: Element {} start must be a single positive integer.",
                label,
                i + 1
            ));
        }
        Some((s[0] as usize) - 1) // 0-based
    };

    let end_usize = if end.is_null() {
        None
    } else {
        let e = end
            .as_integer_slice()
            .ok_or_else(|| anyhow!("{}: Element {} end must be an integer.", label, i + 1))?;
        if e.len() != 1 || e[0] <= 0 {
            return Err(anyhow!(
                "{}: Element {} end must be a single positive integer.",
                label,
                i + 1
            ));
        }
        Some(e[0] as usize) // 1-based exclusive
    };
    SeqRange::new_checked(start_usize, end_usize)
        .map_err(|e| anyhow!("{}: Element {} {}", label, i + 1, e))
}
