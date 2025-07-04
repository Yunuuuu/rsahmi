use std::io::Write;

use super::fastq::FastqRecord;

#[derive(Debug)]
pub(crate) struct FastaRecord<T> {
    pub(crate) id: T,
    pub(crate) desc: Option<T>,
    pub(crate) seq: T,
}

impl<T> From<FastqRecord<T>> for FastaRecord<T> {
    fn from(value: FastqRecord<T>) -> Self {
        Self {
            id: value.id,
            desc: value.desc,
            seq: value.seq,
        }
    }
}

impl<T> FastaRecord<T> {
    pub(crate) fn new(id: T, desc: Option<T>, seq: T) -> Self {
        Self { id, desc, seq }
    }
}

impl<T: AsRef<[u8]>> FastaRecord<T> {
    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
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

    pub(crate) fn as_ref(&self) -> FastaRecord<&[u8]> {
        FastaRecord {
            id: self.id.as_ref(),
            desc: self.desc.as_ref().map(|d| d.as_ref()),
            seq: self.seq.as_ref(),
        }
    }
}

/// Wrapper struct to include UMI and Barcode in the FASTA Record ID
#[derive(Debug)]
pub struct FastaRecordWithUMIBarcode<T> {
    record: FastaRecord<T>,
    umi: Vec<u8>,
    barcode: Vec<u8>,
}

impl<T> FastaRecordWithUMIBarcode<T> {
    pub fn new(record: FastaRecord<T>, umi: Vec<u8>, barcode: Vec<u8>) -> Self {
        Self {
            record,
            umi,
            barcode,
        }
    }
}

impl<T: AsRef<[u8]>> FastaRecordWithUMIBarcode<T> {
    pub fn write<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        let id = self.record.id.as_ref();
        let desc = self.record.desc.as_ref().map(|d| d.as_ref());
        let seq = self.record.seq.as_ref();
        let mut buffer = Vec::with_capacity(
            id.len()
                + b"@RSAHMI:UMI:".len() // Constant length for the UMI part
                + self.umi.len() // UMI length
                + b":BARCODE:".len() // Constant length for the Barcode part
                + self.barcode.len() // Barcode length
                + b":RSAHMI@".len() // Constant RSAHMI end part
                + desc.map(|d| d.len() + 1).unwrap_or(0) // ' '
                + seq.len()
                + 3, // '>' and 2 * '\n'
        );
        buffer.extend_from_slice(b">");
        buffer.extend_from_slice(id);

        // Add UMI and Barcode information
        buffer.extend_from_slice(b"@RSAHMI:UMI:");
        buffer.extend_from_slice(&self.umi);
        buffer.extend_from_slice(b":BARCODE:");
        buffer.extend_from_slice(&self.barcode);
        buffer.extend_from_slice(b":RSAHMI@");

        // Optionally append description
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

#[cfg(test)]
mod tests {
    use super::*;
    // Test the `write` method of `FastaRecord`
    #[test]
    fn test_fasta_record_write() {
        let record = FastaRecord::new(
            b"seq1".to_vec(),
            Some(b"description".to_vec()),
            b"AGCTTAGCTA".to_vec(),
        );

        // Create a mutable vector to capture the written output
        let mut output = Vec::new();

        // Call the write method
        record
            .write(&mut output)
            .expect("Failed to write FASTA record");

        // Expected output
        let expected_output = b">seq1 description\nAGCTTAGCTA\n";

        // Assert that the output matches the expected output
        assert_eq!(output, expected_output);
    }

    // Test the `write` method of `FastaRecordWithUMIBarcode`
    #[test]
    fn test_fasta_record_with_umi_barcode_write() {
        // Example data for the test
        let record = FastaRecord::new(
            b"seq1".to_vec(),
            Some(b"description".to_vec()),
            b"AGCTTAGCTA".to_vec(),
        );
        let umi = b"UMI1234".to_vec();
        let barcode = b"BC5678".to_vec();

        let wrapped_record = FastaRecordWithUMIBarcode::new(record, umi, barcode);

        // Create a mutable vector to capture the written output
        let mut output = Vec::new();

        // Call the write method
        wrapped_record
            .write(&mut output)
            .expect("Failed to write FASTA record");

        // Expected output (including UMI and Barcode in the header)
        let expected_output =
            b">seq1@RSAHMI:UMI:UMI1234:BARCODE:BC5678:RSAHMI@ description\nAGCTTAGCTA\n";

        // Assert that the output matches the expected output
        assert_eq!(output, expected_output);
    }
}
