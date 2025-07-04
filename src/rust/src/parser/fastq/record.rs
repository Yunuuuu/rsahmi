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
        let id = self.id.as_ref();
        let desc = self.desc.as_ref().map(|d| d.as_ref());
        let seq = self.seq.as_ref();
        let sep = self.sep.as_ref();
        let qual = self.qual.as_ref();
        let mut buffer = Vec::with_capacity(
            id.len()
                + desc.map(|d| d.len() + 1).unwrap_or(0) // ' '
                + seq.len()
                + sep.len()
                + qual.len()
                + 5, // '@' and 4 * '\n'
        );
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

        writer.write_all(&buffer)
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
