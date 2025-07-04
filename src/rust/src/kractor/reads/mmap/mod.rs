mod paired;
mod single;

pub(crate) use paired::mmap_kractor_paired_read;
pub(crate) use single::mmap_kractor_single_read;
