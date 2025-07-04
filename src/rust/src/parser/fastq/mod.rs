mod container;
mod error;
mod reader;
mod record;

pub(crate) use container::FastqContainer;
pub(crate) use error::FastqParseError;
pub(crate) use reader::bytes::{FastqBytesChunkPairedReader, FastqBytesChunkReader};
pub(crate) use reader::slice::{FastqSliceChunkPairedReader, FastqSliceChunkReader};
pub(crate) use record::FastqRecord;
