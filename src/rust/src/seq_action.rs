use anyhow::{anyhow, Error, Result};
use bytes::{BufMut, Bytes, BytesMut};
use extendr_api::prelude::*;
use rustc_hash::FxHashMap as HashMap;

use crate::parser::fastq::FastqRecord;
use crate::seq_range::{check_overlap, SeqRange, SeqRanges};
use crate::seq_tag::*;
use crate::utils::*;

pub(crate) struct SubseqEmbedActions {
    tags: TagRanges,
}

impl SubseqEmbedActions {
    fn new(tags: TagRanges) -> Self {
        Self { tags }
    }

    fn has_action(&self) -> bool {
        self.tags.len() > 0
    }

    fn embed(&self, record: &mut FastqRecord<Bytes>) -> Result<()> {
        if self.has_action() {
            let tag_map = self.tags.map_sequences(&record.seq)?;
            record.desc = Some(make_description(
                &tag_map,
                &record.desc.as_ref().map(|d| d.as_ref()),
            ));
        }
        Ok(())
    }
}

struct SubseqTrimActions {
    ranges: SeqRanges,
}

impl SubseqTrimActions {
    fn new(ranges: SeqRanges) -> Self {
        Self { ranges }
    }

    fn has_action(&self) -> bool {
        self.ranges.len() > 0
    }

    /// Trim the sequence and quality string.
    /// If any range is out of bounds, it will return an error.
    fn trim(&self, record: &mut FastqRecord<Bytes>) -> Result<()> {
        if !self.has_action() {
            return Ok(());
        }
        let len = record.seq.len();
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

        let total_seq_len = keep.iter().map(|(start, end)| end - start).sum();
        let mut trimmed_seq = BytesMut::with_capacity(total_seq_len);
        let mut trimmed_qual = BytesMut::with_capacity(total_seq_len);

        for (start, end) in keep {
            trimmed_seq.extend_from_slice(&record.seq[start .. end]);
            trimmed_qual.extend_from_slice(&record.qual[start .. end]);
        }
        record.seq = trimmed_seq.freeze();
        record.qual = trimmed_qual.freeze();
        Ok(())
    }
}

pub(crate) struct SubseqActions {
    embed: SubseqEmbedActions,
    trim: SubseqTrimActions,
}

impl SubseqActions {
    #[allow(dead_code)]
    pub(crate) fn builder() -> SubseqActionsBuilder {
        SubseqActionsBuilder {
            embed_list: Vec::new(),
            trim_list: Vec::new(),
        }
    }

    pub(crate) fn transform_fastq(&self, record: &mut FastqRecord<Bytes>) -> Result<()> {
        self.embed.embed(record)?;
        self.trim.trim(record)?;

        Ok(())
    }
}

/// Paired-end FASTQ transformation logic using optional tag embedding and trimming.
/// Each read end (read1 or read2) can have its own independent action configuration.
pub(crate) struct SubseqPairedActions {
    actions1: Option<SubseqActions>,
    actions2: Option<SubseqActions>,
}

impl SubseqPairedActions {
    pub(crate) fn new(actions1: Option<SubseqActions>, actions2: Option<SubseqActions>) -> Self {
        Self { actions1, actions2 }
    }

    pub(crate) fn transform_fastq(
        &self,
        record1: &mut FastqRecord<Bytes>,
        record2: &mut FastqRecord<Bytes>,
    ) -> Result<()> {
        // Apply tag embedding into `desc` fields
        self.embedded_labels(record1, record2)?;

        // Apply trimming logic for read1
        if let Some(actions) = &self.actions1 {
            actions.trim.trim(record1)?;
        }

        // Apply trimming logic for read2
        if let Some(actions) = &self.actions2 {
            actions.trim.trim(record2)?;
        }
        Ok(())
    }

    /// Embeds tag labels extracted from read1, read2, or both into FASTQ description fields.
    ///
    /// Behavior:
    /// - If only one side has embedding, only that side contributes tags.
    /// - If both sides have embedding, tags are merged by key. If the same tag exists in both,
    ///   the sequences are concatenated in read1-first, read2-second order.
    ///
    /// Tags are serialized using `make_description()` and applied to both reads' description fields.
    fn embedded_labels(
        &self,
        record1: &mut FastqRecord<Bytes>,
        record2: &mut FastqRecord<Bytes>,
    ) -> Result<()> {
        let tag_map = match (&self.actions1, &self.actions2) {
            (Some(actions), None) => actions.embed.tags.map_sequences(&record1.seq)?,
            (None, Some(actions)) => actions.embed.tags.map_sequences(&record2.seq)?,
            (Some(actions1), Some(actions2)) => {
                let mut tag_map = actions1.embed.tags.map_sequences(&record1.seq)?;
                let tag_map2 = actions2.embed.tags.map_sequences(&record2.seq)?;

                // Merge tag→sequence entries
                for (tag, sequence) in tag_map2 {
                    if let Some(v) = tag_map.get_mut(&tag) {
                        v.extend(sequence); // read1 first, read2 second
                    } else {
                        tag_map.insert(tag, sequence);
                    }
                }
                tag_map
            }
            (None, None) => return Ok(()),
        };

        // Only write to description fields if any tag was collected
        if tag_map.len() > 0 {
            record1.desc = Some(make_description(
                &tag_map,
                &record1.desc.as_ref().map(|d| d.as_ref()),
            ));
            record2.desc = Some(make_description(
                &tag_map,
                &record2.desc.as_ref().map(|d| d.as_ref()),
            ));
        }
        Ok(())
    }
}

/// Builder pattern for constructing `SubseqActions` step-by-step.
/// Allows the user to accumulate multiple `embed` and `trim` operations,
/// validating them before finalizing into a concrete `SubseqActions` instance.
pub(crate) struct SubseqActionsBuilder {
    embed_list: Vec<(Bytes, SeqRanges)>,
    trim_list: Vec<SeqRanges>,
}

impl SubseqActionsBuilder {
    fn new() -> Self {
        Self {
            embed_list: Vec::new(),
            trim_list: Vec::new(),
        }
    }

    /// Adds a new embedded tag action with range validation.
    /// Ensures that ranges are sorted and non-overlapping before storing.
    fn add_embed(&mut self, tag: Bytes, mut ranges: SeqRanges) -> Result<()> {
        ranges.sort();
        check_overlap(&ranges)?; // Ensure no conflicting tag extraction ranges
        self.embed_list.push((tag, ranges));
        Ok(())
    }

    /// Adds a trimming action range set. No validation is done here;
    /// validation is deferred to `build()`.
    fn add_trim(&mut self, ranges: SeqRanges) {
        self.trim_list.push(ranges);
    }

    /// Adds a compound action to the builder.
    /// Dispatches to `embed`, `trim`, or both depending on action variant.
    pub(crate) fn add_action(&mut self, action: SeqAction, ranges: SeqRanges) -> Result<()> {
        match action {
            SeqAction::Embed(tag) => {
                self.add_embed(tag, ranges)?;
            }
            SeqAction::Trim => {
                self.add_trim(ranges);
            }
            SeqAction::EmbedTrim(tag) => {
                self.add_embed(tag, ranges.clone())?;
                self.add_trim(ranges);
            }
        }
        Ok(())
    }

    /// Finalizes and builds a `SubseqActions` instance from the accumulated steps.
    /// Sorting and validation of trimming ranges happens here.
    /// Why we sort here:
    /// Trimming ranges from multiple calls may arrive unsorted and overlapping across entries.
    /// Sorting at build time ensures that the full set of ranges is globally ordered and
    /// ready for efficient trimming operations.
    pub(crate) fn build(self) -> Result<SubseqActions> {
        let tags = self.embed_list.into_iter().collect::<TagRanges>();
        let embed_actions = SubseqEmbedActions::new(tags);

        // Flatten all trim ranges into a single sorted list
        let mut full_ranges: SeqRanges = self.trim_list.into_iter().flatten().collect();
        full_ranges.sort();

        // Construct the final SubseqActions struct
        Ok(SubseqActions {
            embed: embed_actions,
            trim: SubseqTrimActions::new(full_ranges),
        })
    }
}

pub(crate) enum SeqAction {
    Embed(Bytes),
    Trim,
    EmbedTrim(Bytes),
}

// Create object from R
pub(crate) fn robj_to_seq_actions<'r>(ranges: &Robj) -> Result<Option<SubseqActions>> {
    if ranges.is_null() {
        return Ok(None);
    }
    Ok(Some(SubseqActions::try_from(ranges)?))
}

impl TryFrom<&Robj> for SubseqActions {
    type Error = Error;
    fn try_from(value: &Robj) -> Result<Self> {
        let mut out = SubseqActionsBuilder::new();
        for ref robj in value
            .as_list()
            .ok_or(anyhow!("Expected a list of sequence range objects."))?
            .values()
        {
            let ranges = SeqRanges::try_from(robj)?;
            let action = SeqAction::try_from(robj)?;
            out.add_action(action, ranges)?;
        }
        out.build()
    }
}

impl TryFrom<&Robj> for SeqAction {
    type Error = Error;
    fn try_from(value: &Robj) -> Result<Self> {
        let resolved_action = if value.inherits("rsahmi_trim") {
            SeqAction::Trim
        } else if value.inherits("rsahmi_embed") {
            SeqAction::Embed(extract_tag_name(value)?)
        } else if value.inherits("rsahmi_embed_trim") {
            SeqAction::EmbedTrim(extract_tag_name(value)?)
        } else {
            return Err(anyhow!(
                "The object does not inherit a valid action class (expected one of: 'rsahmi_trim', 'rsahmi_embed', or 'rsahmi_embed_trim')."
            ));
        };
        Ok(resolved_action)
    }
}

fn make_description(tag_map: &HashMap<Bytes, Vec<&[u8]>>, desc: &Option<&[u8]>) -> Bytes {
    // add prefix, tag and seprator
    let suffix = b'}';
    let mut out = BytesMut::with_capacity(
        // original description length
        desc.map_or(0, |d| d.len() + 1)
            // prefix
            + TAG_PREFIX.len()
            // all tag
            + tag_map.iter().map(|(tag, sequence)| tag.len() + 1 + sequence.iter().map(|s| s.len()).sum::<usize>()).sum::<usize>()
            // all seprator
            + tag_map.len().saturating_sub(1)
            // suffix
            + 1,
    );
    if let Some(v) = desc {
        out.extend_from_slice(&v);
        out.put_u8(b' ');
    }
    out.extend_from_slice(TAG_PREFIX);
    for (i, (tag, sequences)) in tag_map.iter().enumerate() {
        if i > 0 {
            out.put_u8(b':');
        }
        out.extend_from_slice(&tag);
        out.put_u8(b':');
        for seq in sequences {
            out.extend_from_slice(seq);
        }
    }
    out.put_u8(suffix);
    out.freeze()
}
