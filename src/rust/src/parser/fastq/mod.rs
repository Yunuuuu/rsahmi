mod container;
mod error;
mod reader;
mod record;

pub use container::FastqContainer;
pub use error::FastqParseError;
pub use reader::{FastqPairedReader, FastqReader, FastqSource};
#[allow(unused_imports)]
pub use record::FastqRecord;
