mod paired;
mod reader;
mod single;
mod ubread;

pub use paired::mmap_kractor_paired_read;
pub use single::mmap_kractor_single_read;
pub use ubread::mmap_kractor_ubread_read;
