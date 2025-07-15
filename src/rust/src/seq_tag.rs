use anyhow::{anyhow, Error, Result};
use bytes::Bytes;
use extendr_api::prelude::*;
use rustc_hash::FxHashMap as HashMap;

use crate::seq_range::{check_overlap, SeqRanges};

/// A collection of (tag name → sequence ranges) mappings.
/// Each tag (as `Bytes`) maps to a `SeqRanges` defining subsequence locations to extract.
/// Internally stored in a HashMap, and tag order is not guaranteed unless you use an ordered map.
///
/// ## ⚠️ Important:
/// - Each `SeqRanges` **must be sorted and non-overlapping**.
/// - This is **automatically enforced** in the `TryFrom<Robj>` constructor (used in R bindings).
/// - If constructed manually via `TagRanges::new()`, it is the **caller’s responsibility**
///   to ensure `SeqRanges` are **sorted** (`ranges.sort()`) and **non-overlapping** (`check_overlap(&ranges)`).
pub(crate) struct TagRanges {
    map: HashMap<Bytes, SeqRanges>,
}

/// Extracts a tag name from an R object (Robj) used in action annotation.
/// This is used in R interface bindings for embedding.
///
/// # Errors
/// Returns an error if the `"tag"` attribute is missing or not a string.
pub(crate) fn extract_tag_name(robj: &Robj) -> Result<Bytes> {
    robj.get_attrib("tag")
        .and_then(|t| t.as_str())
        .ok_or(anyhow!("'tag' attribute must be provided"))
        .map(|t| Bytes::copy_from_slice(t.as_bytes()))
}

impl TagRanges {
    /// Creates a new `TagRanges` from an existing `HashMap`.
    ///
    /// ⚠️ **Important**: The caller must ensure that all `SeqRanges`:
    /// - Are **sorted** (use `.sort()`).
    /// - Are **non-overlapping** (use `check_overlap(&ranges)`).
    ///
    /// This method does **not** enforce the invariants.
    pub(crate) fn new(tag_ranges: HashMap<Bytes, SeqRanges>) -> Self {
        Self { map: tag_ranges }
    }

    pub(crate) fn len(&self) -> usize {
        self.map.len()
    }

    /// Extract tagged subsequences as slices (`&[u8]`) for each tag.
    ///
    /// This method does not perform any allocation or copying of sequence data.
    /// It simply collects references (`&[u8]`) to the original `seq` slice based
    /// on the defined subsequence ranges for each tag.
    ///
    /// This is useful when:
    /// - You want to inspect or process subsequences without allocating.
    /// - You plan to concatenate or manipulate sequences later (e.g., into `BytesMut`).
    ///
    /// # Returns
    /// - A `HashMap<Bytes, Vec<&[u8]>>`, where each entry represents a tag and its
    ///   corresponding list of extracted subsequence slices from `seq`.
    ///
    /// # Errors
    /// - If any defined `SeqRange` is out of bounds with respect to `seq`, a range error is returned.
    #[inline]
    pub(crate) fn map_sequences<'seq>(
        &self,
        seq: &'seq [u8],
    ) -> Result<HashMap<Bytes, Vec<&'seq [u8]>>> {
        self.map
            .iter()
            .map(|(tag, ranges)| -> Result<(Bytes, Vec<&'seq [u8]>)> {
                // Attempt to extract each defined range from `seq`, collecting borrowed slices
                let sequence = ranges
                    .iter()
                    .map(|r| r.try_extract(seq)) // may return an error
                    .collect::<Result<Vec<_>>>()?;

                // Associate the extracted sequence list with the tag
                Ok((tag.clone(), sequence))
            })
            .collect::<Result<HashMap<Bytes, Vec<&'seq [u8]>>>>()
    }
}

impl IntoIterator for TagRanges {
    type IntoIter = std::collections::hash_map::IntoIter<Bytes, SeqRanges>;
    type Item = (Bytes, SeqRanges);

    fn into_iter(self) -> Self::IntoIter {
        self.map.into_iter()
    }
}

impl<'a> IntoIterator for &'a TagRanges {
    type IntoIter = std::collections::hash_map::Iter<'a, Bytes, SeqRanges>;
    type Item = (&'a Bytes, &'a SeqRanges);

    fn into_iter(self) -> Self::IntoIter {
        (&self.map).into_iter()
    }
}

// Create object from R
pub(crate) fn robj_to_tag_ranges<'r>(ranges: &Robj) -> Result<Option<TagRanges>> {
    if ranges.is_null() {
        return Ok(None);
    }
    Ok(Some(TagRanges::try_from(ranges)?))
}

impl TryFrom<&Robj> for TagRanges {
    type Error = Error;
    fn try_from(value: &Robj) -> Result<Self> {
        let list = value
            .as_list()
            .ok_or(anyhow!("Expected a list of sequence range objects."))?;
        let tag_ranges = list
            .values()
            .into_iter()
            .map(|robj| -> Result<(Bytes, SeqRanges)> {
                if !robj.inherits("rsahmi_tag") {
                    return Err(anyhow!(
                        "The object does not inherit a valid tag class (expected 'rsahmi_tag')."
                    ));
                }
                let tag = extract_tag_name(&robj)?;
                let mut ranges = SeqRanges::try_from(&robj)?;
                ranges.sort();
                check_overlap(&ranges)?;
                Ok((tag, ranges))
            })
            .collect::<Result<HashMap<Bytes, SeqRanges>>>()?;
        Ok(Self::new(tag_ranges))
    }
}
