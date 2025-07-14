use std::io::Write;

#[derive(Debug)]
pub(crate) struct FastqRecord<T> {
    pub(crate) id: T,
    pub(crate) desc: Option<T>,
    pub(crate) seq: T,
    pub(crate) sep: T,
    pub(crate) qual: T,
}

impl<T> FastqRecord<T> {
    #[allow(dead_code)]
    pub(crate) fn new(id: T, desc: Option<T>, seq: T, sep: T, qual: T) -> Self {
        Self {
            id,
            desc,
            seq,
            sep,
            qual,
        }
    }
}

impl<T: AsRef<[u8]>> FastqRecord<T> {
    #[allow(dead_code)]
    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.as_vec())
    }

    pub(crate) fn bytes_size(&self) -> usize {
        self.id.as_ref().len()
            // extra one for space between id and description
            + self.desc.as_ref().map(|d| d.as_ref().len() + 1).unwrap_or(0) // ' '
            + self.seq.as_ref().len()
            + self.sep.as_ref().len()
            + self.qual.as_ref().len()
            + 5 // '@' and 4 * '\n'
    }

    /// Efficiently appends the FASTQ record to the provided Vec<u8>
    pub(crate) fn extend(&self, buf: &mut Vec<u8>) {
        buf.push(b'@');
        buf.extend_from_slice(self.id.as_ref());

        if let Some(desc) = &self.desc {
            buf.push(b' ');
            buf.extend_from_slice(desc.as_ref());
        }

        buf.push(b'\n');
        buf.extend_from_slice(self.seq.as_ref());
        buf.push(b'\n');
        buf.extend_from_slice(self.sep.as_ref());
        buf.push(b'\n');
        buf.extend_from_slice(self.qual.as_ref());
        buf.push(b'\n');
    }

    pub(crate) fn as_vec(&self) -> Vec<u8> {
        let id = self.id.as_ref();
        let desc = self.desc.as_ref().map(|d| d.as_ref());
        let seq = self.seq.as_ref();
        let sep = self.sep.as_ref();
        let qual = self.qual.as_ref();
        let mut buffer = Vec::with_capacity(self.bytes_size());
        buffer.push(b'@');
        buffer.extend_from_slice(id);

        if let Some(desc) = &desc {
            buffer.push(b' ');
            buffer.extend_from_slice(desc);
        }

        buffer.push(b'\n');
        buffer.extend_from_slice(seq);
        buffer.push(b'\n');
        buffer.extend_from_slice(sep);
        buffer.push(b'\n');
        buffer.extend_from_slice(qual);
        buffer.push(b'\n');
        buffer
    }

    #[allow(dead_code)]
    pub(crate) fn write_buf(&self, buf: &mut [u8]) -> std::io::Result<usize> {
        let mut pos = 0;

        // Write '@' and ID
        buf[pos] = b'@';
        pos += 1;

        let id = self.id.as_ref();
        if pos + id.len() > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + id.len()].copy_from_slice(id);
        pos += id.len();

        // Optional description
        if let Some(desc) = &self.desc {
            if pos + 1 > buf.len() {
                return Err(std::io::ErrorKind::WriteZero.into());
            }
            buf[pos] = b' ';
            pos += 1;

            let desc = desc.as_ref();
            if pos + desc.len() > buf.len() {
                return Err(std::io::ErrorKind::WriteZero.into());
            }
            buf[pos .. pos + desc.len()].copy_from_slice(desc);
            pos += desc.len();
        }

        // Line break
        if pos + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos] = b'\n';
        pos += 1;

        // Sequence
        let seq = self.seq.as_ref();
        if pos + seq.len() + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + seq.len()].copy_from_slice(seq);
        pos += seq.len();
        buf[pos] = b'\n';
        pos += 1;

        // Separator
        let sep = self.sep.as_ref();
        if pos + sep.len() + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + sep.len()].copy_from_slice(sep);
        pos += sep.len();
        buf[pos] = b'\n';
        pos += 1;

        // Quality
        let qual = self.qual.as_ref();
        if pos + qual.len() + 1 > buf.len() {
            return Err(std::io::ErrorKind::WriteZero.into());
        }
        buf[pos .. pos + qual.len()].copy_from_slice(qual);
        pos += qual.len();
        buf[pos] = b'\n';
        pos += 1;

        Ok(pos)
    }

    #[allow(dead_code)]
    pub(crate) fn as_ref(&self) -> FastqRecord<&[u8]> {
        FastqRecord {
            id: self.id.as_ref(),
            desc: self.desc.as_ref().map(|d| d.as_ref()),
            seq: self.seq.as_ref(),
            sep: self.sep.as_ref(),
            qual: self.qual.as_ref(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*; // Import FastqRecord

    #[test]
    fn test_fastq_record_write_with_description() {
        let record = FastqRecord::new(
            b"SEQ_ID".as_ref(),
            Some(b"desc".as_ref()),
            b"ACGTACGT".as_ref(),
            b"+".as_ref(),
            b"IIIIIIII".as_ref(),
        );

        let mut output = Cursor::new(Vec::new());
        record.write(&mut output).expect("Write failed");

        let expected = b"@SEQ_ID desc\nACGTACGT\n+\nIIIIIIII\n";
        assert_eq!(output.into_inner(), expected);
    }

    #[test]
    fn test_fastq_record_write_without_description() {
        let record = FastqRecord::new(
            b"SEQ_ID".as_ref(),
            None,
            b"ACGTACGT".as_ref(),
            b"+".as_ref(),
            b"IIIIIIII".as_ref(),
        );

        let mut output = Cursor::new(Vec::new());
        record.write(&mut output).expect("Write failed");

        let expected = b"@SEQ_ID\nACGTACGT\n+\nIIIIIIII\n";
        assert_eq!(output.into_inner(), expected);
    }
}
