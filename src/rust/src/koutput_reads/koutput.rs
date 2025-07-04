use std::io::Read;

use aho_corasick::AhoCorasick;
use anyhow::{anyhow, Result};
use bytes::Bytes;
use crossbeam_channel::{Receiver, Sender};
use memchr::{memchr, memchr2};
use rustc_hash::FxHashMap as HashMap;

use crate::batchsender::BatchSender;
use crate::kractor::koutput::io::KoutputBytesChunkReader;
use crate::kractor::koutput::KOUTPUT_TAXID_PREFIX_FINDER;
use crate::reader::bytes::BytesProgressBarReader;

pub(crate) fn reader_parse_koutput<R: Read + Send>(
    reader: BytesProgressBarReader<R>,
    include_aho: AhoCorasick,
    exclude_aho: Option<AhoCorasick>,
    chunk_size: usize,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<HashMap<Bytes, (Bytes, Bytes, Bytes)>> {
    // for kmer, we counts total and unique k-mers per taxon across cell barcodes,
    // using both the cell barcode and unique molecular identifier (UMI) to resolve
    // read identity at the single-cell level. It aggregates k-mer counts for each
    // taxonomic rank of interest (by default, genus and species), including all
    // descendant taxa within those ranks.
    std::thread::scope(|scope| {
        // Create a channel between the parser and writer threads
        // The channel transmits batches
        let (koutput_tx, koutput_rx): (
            Sender<Vec<(Bytes, (Bytes, Bytes, Bytes))>>,
            Receiver<Vec<(Bytes, (Bytes, Bytes, Bytes))>>,
        ) = crate::new_channel(nqueue);

        // ─── Parser Thread ─────────────────────────────────────
        // Streams Kraken2 output data, filters by ID set
        let mut reader = KoutputBytesChunkReader::with_capacity(chunk_size, reader);
        let koutput_handle = scope.spawn(move || -> Result<()> {
            // will move `reader`, `koutput_tx`, and `matcher`
            rayon::scope(|s| -> Result<()> {
                while let Some(chunk) = reader.chunk_reader()? {
                    s.spawn(|_| {
                        let mut chunk = chunk; // Move the chunk
                        let mut thread_tx =
                            BatchSender::with_capacity(batch_size, koutput_tx.clone());
                        'chunk_loop: while let Some(line) = chunk.read_line() {
                            let mut field_start = 0usize;
                            let mut field_index = 0usize;
                            let mut sequence_id = None;
                            let mut taxid = None;
                            let lca;
                            while let Some(tab_pos) = memchr(b'\t', &line[field_start ..]) {
                                let field = &line[field_start .. (field_start + tab_pos)];
                                if field_index == 0 {
                                    // Field 0: "C" or "U" — skip unclassified reads
                                    if field.len() != 1 || field[0] != b'C' {
                                        continue 'chunk_loop;
                                    }
                                } else if field_index == 1 {
                                    // Save sequence_id field (field 2)
                                    sequence_id = Some(field);
                                } else if field_index == 2 {
                                    // Save taxid field (field 3) if it passes filtering
                                    if let Some(start) = KOUTPUT_TAXID_PREFIX_FINDER.find(field) {
                                        let mut input = aho_corasick::Input::new(field);
                                        input.set_start(start);
                                        // Skip this line if taxid is not in matcher
                                        if include_aho.find(field).is_none() {
                                            continue 'chunk_loop;
                                        };
                                        taxid = Some(field);
                                    } else {
                                        // Skip line if taxid doesn't contain the prefix
                                        continue 'chunk_loop;
                                    }
                                } else if field_index == 3 {
                                    // Field 4 (LCA): remainder after the last tab
                                    field_start += tab_pos + 1;
                                    // Ignore final \n character
                                    if let Some(pos) = memchr2(b'\n', b'\t', &line[field_start ..])
                                    {
                                        lca = &line[field_start .. (field_start + pos)]
                                    } else {
                                        lca = &line[field_start ..]
                                    };
                                    if let Some(ref exclude_matcher) = exclude_aho {
                                        if exclude_matcher.find(lca).is_some() {
                                            continue 'chunk_loop;
                                        }
                                    }
                                    if let (Some(sequence_id), Some(taxid)) = (sequence_id, taxid) {
                                        // Just ignore the error message, it won't occur
                                        let _ = thread_tx.send((
                                            line.slice_ref(sequence_id),
                                            (
                                                line.slice_ref(field), // sequence length
                                                line.slice_ref(taxid),
                                                line.slice_ref(lca),
                                            ),
                                        ));
                                    };
                                    continue 'chunk_loop;
                                }
                                field_index += 1;
                                field_start += tab_pos + 1;
                            }
                        }
                        // Just ignore the error message, it won't occur
                        let _ = thread_tx.flush();
                    });
                }
                Ok(())
            })
        });
        koutput_handle
            .join()
            .map_err(|e| anyhow!("Parsing Koutput thread panicked: {:?}", e))??;
        Ok(koutput_rx
            .into_iter()
            .flatten()
            .collect::<HashMap<Bytes, (Bytes, Bytes, Bytes)>>())
    })
}
