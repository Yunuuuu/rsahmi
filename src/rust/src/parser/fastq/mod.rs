mod error;
mod reader;
mod record;
mod container;

pub use error::FastqParseError;
pub use reader::{FastqReader, FastqSource, FastqPairedReader};
#[allow(unused_imports)]
pub use record::FastqRecord;
pub use container::FastqContainer;
