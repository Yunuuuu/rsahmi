mod container;
mod error;
mod record;

pub use container::FastqContainer;
pub use error::FastqParseError;
#[allow(unused_imports)]
pub use record::FastqRecord;
