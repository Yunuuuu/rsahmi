use std::path::Path;

use anyhow::{anyhow, Context, Result};
use bytes::{Bytes, BytesMut};
use crossbeam_channel::{Receiver, Sender};
use indicatif::{ProgressBar, ProgressFinish};
use memchr::memchr;
use memchr::memmem::Finder;
use rustc_hash::FxHashMap as HashMap;
use rustc_hash::FxHashSet as HashSet;

mod counter;

use counter::{CountTotal, CountUnique, Countable};

use crate::batchsender::BatchSender;
use crate::reader::LineReader;
use crate::utils::*;

/// Return `true` if all base counts are ≤ `threshold`, otherwise `false`.
fn pass_complexity_filter(seq: &[u8], threshold: usize) -> bool {
    // remove low complexity reads (<20 non-sequentially repeated nucleotides)
    let threshold = seq.len() - threshold;
    let mut counts = HashMap::with_capacity_and_hasher(4, rustc_hash::FxBuildHasher); // ATGC
    for &b in seq {
        let count = counts.entry(b).or_insert(0);
        *count += 1;
        if *count > threshold {
            return false; // Early exit: too many repeats
        }
    }
    true
}

/// Returns `true` if all quality scores are ≥ `min_phred`.
fn pass_quality_filter(qual: &[u8], threshold: u8) -> bool {
    // threshold 53 for Phred score < 20 (Phred+33 ASCII)
    // threshold 84 for Phred score < 20 (Phred+64 ASCII)
    qual.iter().all(|&q| q >= threshold)
}

/// ReadCounter tracks reads either by total count or unique UMI.
enum ReadCounter {
    /// Count all reads without considering UMI (bulk mode).
    Total(CountTotal),
    /// Count only unique UMIs (single-cell mode).
    Unique(CountUnique<Bytes>),
}

impl ReadCounter {
    fn insert(&mut self, item: Option<&[u8]>) {
        match self {
            ReadCounter::Total(inner) => inner.insert(()),
            ReadCounter::Unique(inner) => {
                // SAFETY: `item` is guaranteed to be `Some` in Unique mode,
                // as enforced by logic in the `sckmer()` function.
                inner.insert(Bytes::copy_from_slice(unsafe { item.unwrap_unchecked() }))
            }
        }
    }

    fn count(&self) -> usize {
        match self {
            ReadCounter::Total(inner) => inner.count(),
            ReadCounter::Unique(inner) => inner.count(),
        }
    }
}

/// ReadsAndKmer holds per-(barcode, taxon) statistics:
/// number of reads, total k-mers, and unique k-mers.
pub(super) struct ReadsAndKmer {
    reads: ReadCounter,
    kmer_total: CountTotal,
    kmer_unique: CountUnique<Bytes>,
}

impl ReadsAndKmer {
    fn new(reads: ReadCounter, kmer_total: CountTotal, kmer_unique: CountUnique<Bytes>) -> Self {
        Self {
            reads,
            kmer_total,
            kmer_unique,
        }
    }

    pub(super) fn reads(&self) -> usize {
        self.reads.count()
    }

    pub(super) fn kmer_total(&self) -> usize {
        self.kmer_total.count()
    }

    pub(super) fn kmer_unique(&self) -> usize {
        self.kmer_unique.count()
    }

    fn add_read(&mut self, umi: Option<&[u8]>) {
        self.reads.insert(umi);
    }

    /// Extract and add k-mers from a sequence and its LCA annotation.
    #[allow(dead_code)]
    fn add_kmer(&mut self, lca: &[u8], sequence: &[u8]) -> Result<()> {
        self.add_kmers(&extract_kmers(lca, sequence)?);
        Ok(())
    }

    /// Insert a batch of k-mers, updating total and unique counts.
    fn add_kmers(&mut self, kmers: &[Bytes]) {
        for kmer in kmers {
            self.kmer_total.insert(());
            self.kmer_unique.insert(kmer.clone());
        }
    }
}

/// Parses a Koutreads-format file and counts reads and k-mers per (barcode, taxon).
/// Each taxon aggregates k-mers from its descendant taxa. Optionally groups reads
/// by barcode and/or UMI if tags are provided.
pub(super) fn count_kmers_and_reads<'taxid, P: AsRef<Path> + ?Sized>(
    koutreads: &P,
    ancestor_map: HashMap<&[u8], HashSet<&'taxid [u8]>>,
    umi_tag: Option<&str>,
    barcode_tag: Option<&str>,
    batch_size: usize,
    nqueue: Option<usize>,
) -> Result<HashMap<Bytes, HashMap<&'taxid [u8], ReadsAndKmer>>> {
    let input: &Path = koutreads.as_ref();
    let style = progress_reader_style()?;
    let pb = ProgressBar::new(input.metadata()?.len() as u64).with_finish(ProgressFinish::Abandon);
    pb.set_prefix("Parsing Koutreads");
    pb.set_style(style);

    // This function processes a Koutreads-format file and collects k-mer counts
    // per (barcode, taxon). Each taxon aggregates k-mers from its descendant taxa.
    // Reads can optionally be grouped by UMI and/or barcode tags.
    std::thread::scope(
        |scope| -> Result<HashMap<Bytes, HashMap<&[u8], ReadsAndKmer>>> {
            // Shared queue between reader and parser threads
            let (reader_tx, reader_rx): (Sender<Vec<BytesMut>>, Receiver<Vec<BytesMut>>) =
                new_channel(nqueue);

            // ─── Parser Thread ─────────────────────────────────────
            // Consumes batches of lines, parses fields, extracts barcode/UMI/LCA/kmers,
            // and accumulates stats into (barcode, taxon) → SCKmer map
            let parser_handle = scope.spawn(
                move || -> Result<HashMap<Bytes, HashMap<&[u8], ReadsAndKmer>>> {
                    let mut barcode_taxon_map =
                        HashMap::with_capacity_and_hasher(1, rustc_hash::FxBuildHasher);
                    let umi_finder = umi_tag.as_ref().map(|tag| Finder::new(tag));
                    let barcode_finder = barcode_tag.as_ref().map(|tag| Finder::new(tag));

                    while let Ok(lines) = reader_rx.recv() {
                        for line in lines {
                            let line = line.freeze();
                            let fields: Vec<&[u8]> = line.split(|b| *b == b'\t').collect();
                            if fields.len() != 5 {
                                return Err(anyhow!("Invalid file: must have 5 fields"));
                            }

                            // ─── Extract and validate fields ───────────────
                            // taxid + tags + lca + seq + qual
                            let qual = unsafe { fields.get_unchecked(4) };
                            if !pass_quality_filter(qual, 53) {
                                continue;
                            }
                            let seq = unsafe { fields.get_unchecked(3) };
                            if !pass_complexity_filter(seq, 20) {
                                continue;
                            }
                            let taxid = unsafe { fields.get_unchecked(0) };

                            // ─── Resolve taxonomic ancestors ───────────────
                            if let Some(ancestors) = ancestor_map.get(taxid) {
                                // ─── Extract barcode and UMI (optional) ────────
                                let tags = unsafe { fields.get_unchecked(1) };
                                let barcode = extract_tag(tags, &barcode_finder, &barcode_tag)
                                    .with_context(|| {
                                        format!(
                                            "Failed to extract barcode in line '{}'",
                                            String::from_utf8_lossy(&line)
                                        )
                                    })?;
                                let umi = extract_tag(tags, &umi_finder, &umi_tag).with_context(
                                    || {
                                        format!(
                                            "Failed to extract umi in line '{}'",
                                            String::from_utf8_lossy(&line)
                                        )
                                    },
                                )?;

                                let barcode = barcode
                                    .map(Bytes::copy_from_slice)
                                    .unwrap_or_else(Bytes::new); // Default: treat as single-cell
                                let barcode_map =
                                    barcode_taxon_map.entry(barcode).or_insert_with(|| {
                                        HashMap::with_capacity_and_hasher(
                                            1,
                                            rustc_hash::FxBuildHasher,
                                        )
                                    });

                                // ─── Extract all kmers from sequence(s) ─────
                                // A space-delimited list indicating the LCA mapping of each
                                // k-mer in the sequence(s). For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
                                //
                                // the first 13 k-mers mapped to taxonomy ID #562
                                // the next 4 k-mers mapped to taxonomy ID #561
                                // the next 31 k-mers contained an ambiguous nucleotide
                                // the next k-mer was not in the database
                                // the last 3 k-mers mapped to taxonomy ID #562
                                let lca = unsafe { fields.get_unchecked(2) };
                                let kmers =
                                    match (LCA_SEPARATOR_FINDER.find(lca), memchr(b' ', seq)) {
                                        (Some(lca_pos), Some(seq_pos)) => {
                                            // Paired-end
                                            // Note that paired read data will contain a "|:|" token in this
                                            // list to indicate the end of one read and the beginning of another.
                                            let lca1 = &lca[.. lca_pos];
                                            let lca2 = &lca[lca_pos + LCA_SEPARATOR.len() + 1 ..];
                                            let seq1 = &seq[.. seq_pos];
                                            let seq2 = &seq[seq_pos + 2 ..];
                                            [extract_kmers(lca1, seq1)?, extract_kmers(lca2, seq2)?]
                                                .concat()
                                        }
                                        (None, None) => {
                                            // Single-end
                                            extract_kmers(lca, seq)?
                                        }
                                        (_, _) => {
                                            return Err(anyhow!("Mismatched LCA/sequence format"));
                                        }
                                    };

                                // ─── Update stats per (barcode, ancestor taxon) ───────
                                for ancestor in ancestors {
                                    let entry = barcode_map.entry(*ancestor).or_insert_with(|| {
                                        let reads = if umi_tag.is_some() {
                                            ReadCounter::Unique(CountUnique::with_capacity(1))
                                        } else {
                                            ReadCounter::Total(CountTotal::new())
                                        };
                                        ReadsAndKmer::new(
                                            reads,
                                            CountTotal::new(),
                                            CountUnique::with_capacity(1),
                                        )
                                    });
                                    entry.add_read(umi);
                                    entry.add_kmers(&kmers);
                                }
                            }
                        }
                    }
                    Ok(barcode_taxon_map)
                },
            );

            // ─── reader Thread ─────────────────────────────────────
            // Reads lines from input file and sends them in batches to parser thread
            let reader_handle = scope.spawn(move || -> Result<()> {
                let mut reader = LineReader::with_capacity(
                    BUFFER_SIZE,
                    new_reader(input, BUFFER_SIZE, Some(pb))?,
                );
                let mut reader_tx: BatchSender<BytesMut> =
                    BatchSender::with_capacity(batch_size, reader_tx);
                while let Some(line) = reader
                    .read_line()
                    .with_context(|| format!("(Reader) Failed to read line"))?
                {
                    if line.iter().all(|b| b.is_ascii_whitespace()) {
                        continue;
                    }
                    reader_tx.send(line).with_context(|| {
                        format!("(Reader) Failed to send lines to Parser thread")
                    })?;
                }
                reader_tx
                    .flush()
                    .with_context(|| format!("(Reader) Failed to flush lines to Parser thread"))?;
                Ok(())
            });

            // ─── Join Threads and Propagate Errors ────────────────
            let out = parser_handle
                .join()
                .map_err(|e| anyhow!("(Parser) thread panicked: {:?}", e))??;
            reader_handle
                .join()
                .map_err(|e| anyhow!("(Reader) thread panicked: {:?}", e))??;
            Ok(out)
        },
    )
}

const LCA_SEPARATOR: &'static [u8] = b"|:|";
static LCA_SEPARATOR_FINDER: std::sync::LazyLock<Finder> =
    std::sync::LazyLock::new(|| Finder::new(TAG_PREFIX));

fn extract_tag<'t>(
    tags: &'t [u8],
    finder: &Option<Finder>,
    label: &Option<&str>,
) -> Result<Option<&'t [u8]>> {
    match (finder, label) {
        (Some(finder), Some(lab)) => {
            if let Some(start) = finder.find(tags) {
                // Skip over the tag label and colon (e.g., "TAG:")
                let start = start + lab.len() + 1;
                if start >= tags.len() {
                    return Err(anyhow!("Tag '{}' found, but no value follows", lab));
                }

                // Look for the next space to determine the end of the tag value
                let end = memchr(b' ', &tags[start ..])
                    .map(|e| start + e)
                    .unwrap_or(tags.len());

                Ok(Some(&tags[start .. end]))
            } else {
                return Err(anyhow!("Tag '{}' not found in input", lab));
            }
        }
        (None, None) => Ok(None),
        _ => unreachable!(),
    }
}

fn extract_kmers(lca: &[u8], sequence: &[u8]) -> Result<Vec<Bytes>> {
    // Step 1: Parse all `taxid:num_kmers` pairs
    let num_kmers = lca
        .trim_ascii()
        .split(|b| *b == b' ')
        .map(|pair| {
            if let Some(pos) = memchr(b':', pair) {
                // SAFETY: checked pos + 1 < pair.len()
                if pos + 1 >= pair.len() {
                    return Err(anyhow!(
                        "Invalid lca pair, missing number after ':' in {:?}",
                        lca
                    ));
                }
                let n = std::str::from_utf8(unsafe { pair.get_unchecked(pos + 1 ..) })?
                    .parse::<usize>()?;
                Ok(n)
            } else {
                Err(anyhow!("Invalid lca pair, missing ':' in {:?}", lca))
            }
        })
        .collect::<Result<Vec<usize>>>()?;

    // Step 2: Compute total k-mers and derive k-mer length
    let total_kmers: usize = num_kmers.iter().sum();
    if total_kmers == 0 || total_kmers > sequence.len() {
        return Err(anyhow!(
            "Invalid total kmer count: {}, sequence length: {}",
            total_kmers,
            sequence.len()
        ));
    }

    let k = sequence.len() + 1 - total_kmers;

    // Step 3: Extract and insert k-mers
    let sequence = Bytes::copy_from_slice(sequence);
    let kmers = sequence.windows(k).map(|w| sequence.slice_ref(w)).collect();

    Ok(kmers)
}
