mod paired;
mod reader;
mod single;
mod ubread;

pub use paired::reader_kractor_paired_read;
pub use single::reader_kractor_single_read;
pub use ubread::reader_kractor_ubread_read;
