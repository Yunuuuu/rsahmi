use std::io::Write;

#[derive(Debug)]
pub struct FastaRecord<T> {
    pub id: T,
    pub desc: Option<T>,
    pub seq: T,
}

impl<T> FastaRecord<T> {
    pub fn new(id: T, desc: Option<T>, seq: T) -> Self {
        Self { id, desc, seq }
    }
}

impl<T: AsRef<[u8]>> FastaRecord<T> {
    pub fn write<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let id = self.id.as_ref();
        let desc = self.desc.as_ref().map(|d| d.as_ref());
        let seq = self.seq.as_ref();
        let mut buffer = Vec::with_capacity(
            id.len()
                + desc.map(|d| d.len() + 1).unwrap_or(0) // ' '
                + seq.len()
                + 3, // '>' and 2 * '\n'
        );
        buffer.extend_from_slice(b">");
        buffer.extend_from_slice(id);

        if let Some(desc) = &desc {
            buffer.push(b' ');
            buffer.extend_from_slice(desc);
        }

        buffer.push(b'\n');
        buffer.extend_from_slice(seq);
        buffer.push(b'\n');

        writer.write_all(&buffer)
    }
}
