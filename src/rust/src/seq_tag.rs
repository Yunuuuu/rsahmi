use anyhow::{anyhow, Result};
use bytes::Bytes;
use extendr_api::prelude::*;
use rustc_hash::FxHashMap as HashMap;

use crate::seq_range::SeqRanges;

/// A collection of (tag name → sequence ranges) mappings.
/// Each tag (as `Bytes`) maps to a `SeqRanges` defining subsequence locations to extract.
/// Internally stored in a HashMap, and tag order is not guaranteed unless you use an ordered map.
pub(crate) struct TagRanges {
    map: HashMap<Bytes, SeqRanges>,
}

/// Extracts a tag name from an R object (Robj) used in action annotation.
/// This is used in R interface bindings for embedding.
///
/// # Errors
/// Returns an error if the `"tag"` attribute is missing or not a string.
pub(crate) fn extract_tag_name(action: &Robj) -> Result<Bytes> {
    action
        .get_attrib("tag")
        .and_then(|t| t.as_str())
        .ok_or(anyhow!("'tag' attribute must be provided"))
        .map(|t| Bytes::copy_from_slice(t.as_bytes()))
}

impl TagRanges {
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

impl FromIterator<(Bytes, SeqRanges)> for TagRanges {
    fn from_iter<T: IntoIterator<Item = (Bytes, SeqRanges)>>(iter: T) -> Self {
        Self {
            map: iter.into_iter().collect::<HashMap<Bytes, SeqRanges>>(),
        }
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
