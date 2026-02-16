//! K-mer counting and output.
//!
//! This module provides the main k-mer counting functionality, using parallel
//! processing to efficiently count canonical k-mers across all sequences in FASTA or FASTQ files.

use crate::{
    cli::OutputFormat,
    format::SequenceFormat,
    input::Input,
    kmer::{pack_canonical, unpack_to_string, KmerLength},
    progress::{Progress, ProgressTracker},
    reader::{read, read_with_quality, SequenceWithQuality},
    streaming::count_kmers_stdin_with_format,
};
use bytes::Bytes;
use dashmap::DashMap;
use rayon::prelude::{ParallelBridge, ParallelIterator};
use rustc_hash::FxHasher;
use serde::Serialize;
use std::{
    collections::HashMap,
    error::Error,
    fmt::Debug,
    hash::BuildHasherDefault,
    io::{stdout, BufWriter, Error as IoError, Write},
    path::Path,
};
use thiserror::Error;

#[cfg(feature = "tracing")]
#[allow(unused_imports)]
use tracing::{debug, info, info_span};

/// Errors that can occur during k-mer processing.
#[derive(Debug, Error)]
pub enum ProcessError {
    /// Error reading input FASTA file.
    #[error("Unable to read input: {0}")]
    ReadError(#[from] Box<dyn Error>),

    /// Error writing output.
    #[error("Unable to write output: {0}")]
    WriteError(#[from] IoError),

    /// Error serializing JSON output.
    #[error("Unable to serialize JSON: {0}")]
    JsonError(#[from] serde_json::Error),
}

/// A k-mer with its count, used for JSON serialization.
#[derive(Serialize)]
struct KmerCount {
    kmer: String,
    count: u64,
}

/// Counts k-mers in a FASTA file and writes results to stdout.
///
/// Reads sequences from the FASTA file at `path`, counts all canonical k-mers
/// of length `k`, and writes the results to stdout in FASTA-like format.
///
/// # Errors
///
/// Returns `ProcessError::ReadError` if the file cannot be read.
/// Returns `ProcessError::WriteError` if output cannot be written.
pub fn run<P>(path: P, k: usize) -> Result<(), ProcessError>
where
    P: AsRef<Path> + Debug,
{
    run_with_options(path, k, OutputFormat::Fasta, 1)
}

/// Counts k-mers with configurable output format and filtering.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
/// * `k` - K-mer length
/// * `format` - Output format (fasta, tsv, or json)
/// * `min_count` - Minimum count threshold (k-mers below this are excluded)
///
/// # Errors
///
/// Returns `ProcessError` on read, write, or serialization errors.
pub fn run_with_options<P>(
    path: P,
    k: usize,
    format: OutputFormat,
    min_count: u64,
) -> Result<(), ProcessError>
where
    P: AsRef<Path> + Debug,
{
    let counts = count_kmers(&path, k)?;
    output_counts(counts, format, min_count)
}

/// Counts k-mers from an input source (file or stdin) and writes results to stdout.
///
/// This is the main entry point for input-agnostic k-mer counting with output.
/// Input format is auto-detected from file extension.
///
/// # Arguments
///
/// * `input` - The input source (file path or stdin)
/// * `k` - K-mer length
/// * `output_format` - Output format (fasta, tsv, or json)
/// * `min_count` - Minimum count threshold (k-mers below this are excluded)
///
/// # Errors
///
/// Returns `ProcessError` on read, write, or serialization errors.
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::run::run_with_input;
/// use kmerust::input::Input;
/// use kmerust::cli::OutputFormat;
/// use std::path::Path;
///
/// // From file
/// let input = Input::from_path(Path::new("genome.fa"));
/// run_with_input(&input, 21, OutputFormat::Tsv, 1)?;
///
/// // From stdin
/// let input = Input::Stdin;
/// // run_with_input(&input, 21, OutputFormat::Tsv, 1)?;
/// # Ok::<(), kmerust::run::ProcessError>(())
/// ```
pub fn run_with_input(
    input: &Input,
    k: usize,
    output_format: OutputFormat,
    min_count: u64,
) -> Result<(), ProcessError> {
    run_with_input_format(input, k, output_format, min_count, SequenceFormat::Auto)
}

/// Counts k-mers from an input source with explicit input format specification.
///
/// # Arguments
///
/// * `input` - The input source (file path or stdin)
/// * `k` - K-mer length
/// * `output_format` - Output format (fasta, tsv, or json)
/// * `min_count` - Minimum count threshold (k-mers below this are excluded)
/// * `input_format` - Input file format (Auto, Fasta, or Fastq)
///
/// # Errors
///
/// Returns `ProcessError` on read, write, or serialization errors.
pub fn run_with_input_format(
    input: &Input,
    k: usize,
    output_format: OutputFormat,
    min_count: u64,
    input_format: SequenceFormat,
) -> Result<(), ProcessError> {
    let counts = match input {
        Input::File(path) => count_kmers_with_format(path, k, input_format)?,
        Input::Stdin => count_kmers_stdin_with_format(k, input_format)
            .map_err(|e| ProcessError::ReadError(e.into()))?,
    };
    output_counts(counts, output_format, min_count)
}

/// Counts k-mers from an input source with quality filtering and writes to stdout.
///
/// For FASTQ input, k-mers containing bases with quality below `min_quality`
/// are skipped. For FASTA input, the quality parameter is ignored.
///
/// # Arguments
///
/// * `input` - The input source (file path or stdin)
/// * `k` - K-mer length
/// * `output_format` - Output format (fasta, tsv, or json)
/// * `min_count` - Minimum count threshold
/// * `input_format` - Input file format
/// * `min_quality` - Optional minimum Phred quality score (0-93)
///
/// # Errors
///
/// Returns `ProcessError` on read, write, or serialization errors.
pub fn run_with_quality(
    input: &Input,
    k: usize,
    output_format: OutputFormat,
    min_count: u64,
    input_format: SequenceFormat,
    min_quality: Option<u8>,
) -> Result<(), ProcessError> {
    let counts = match input {
        Input::File(path) => count_kmers_with_quality(path, k, input_format, min_quality)?,
        // For stdin, quality filtering is not yet supported
        Input::Stdin => count_kmers_stdin_with_format(k, input_format)
            .map_err(|e| ProcessError::ReadError(e.into()))?,
    };
    output_counts(counts, output_format, min_count)
}

/// Counts k-mers and returns them as a `HashMap`.
///
/// This is the main library API for counting k-mers without writing to stdout.
/// Input format is auto-detected from the file extension.
///
/// # Arguments
///
/// * `path` - Path to the FASTA or FASTQ file
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A `HashMap` mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be read
pub fn count_kmers<P>(path: P, k: usize) -> Result<HashMap<String, u64>, Box<dyn Error>>
where
    P: AsRef<Path> + Debug,
{
    count_kmers_with_format(path, k, SequenceFormat::Auto)
}

/// Counts k-mers with explicit format specification.
///
/// # Arguments
///
/// * `path` - Path to the sequence file
/// * `k` - K-mer length (must be 1-32)
/// * `format` - Input file format (Auto, Fasta, or Fastq)
///
/// # Returns
///
/// A `HashMap` mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be read
pub fn count_kmers_with_format<P>(
    path: P,
    k: usize,
    format: SequenceFormat,
) -> Result<HashMap<String, u64>, Box<dyn Error>>
where
    P: AsRef<Path> + Debug,
{
    #[cfg(feature = "tracing")]
    info!(k = k, path = ?path, "Starting k-mer counting");

    // Validate k-mer length upfront to provide a clear error
    let k_len = KmerLength::new(k)?;

    #[cfg(feature = "tracing")]
    let read_span = info_span!("read_sequences", path = ?path).entered();

    let sequences = read(&path, format)?;

    #[cfg(feature = "tracing")]
    drop(read_span);

    #[cfg(feature = "tracing")]
    let process_span = info_span!("process_sequences").entered();

    let kmer_map = KmerMap::new().build(sequences, k);

    #[cfg(feature = "tracing")]
    drop(process_span);

    let result = kmer_map.into_hashmap(k_len);

    #[cfg(feature = "tracing")]
    info!(unique_kmers = result.len(), "K-mer counting complete");

    Ok(result)
}

/// Counts k-mers with quality filtering.
///
/// For FASTQ input, k-mers containing bases with quality below `min_quality`
/// are skipped. For FASTA input, the quality parameter is ignored.
///
/// # Arguments
///
/// * `path` - Path to the FASTA or FASTQ file
/// * `k` - K-mer length (must be 1-32)
/// * `format` - Input file format (Auto, Fasta, or Fastq)
/// * `min_quality` - Optional minimum Phred quality score (0-93)
///
/// # Returns
///
/// A `HashMap` mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be read
pub fn count_kmers_with_quality<P>(
    path: P,
    k: usize,
    format: SequenceFormat,
    min_quality: Option<u8>,
) -> Result<HashMap<String, u64>, Box<dyn Error>>
where
    P: AsRef<Path> + Debug,
{
    #[cfg(feature = "tracing")]
    info!(k = k, path = ?path, min_quality = ?min_quality, "Starting k-mer counting with quality filter");

    // Validate k-mer length upfront to provide a clear error
    let k_len = KmerLength::new(k)?;

    #[cfg(feature = "tracing")]
    let read_span = info_span!("read_sequences", path = ?path).entered();

    let sequences = read_with_quality(&path, format)?;

    #[cfg(feature = "tracing")]
    drop(read_span);

    #[cfg(feature = "tracing")]
    let process_span = info_span!("process_sequences").entered();

    let kmer_map = KmerMap::new().build_with_quality(sequences, k, min_quality);

    #[cfg(feature = "tracing")]
    drop(process_span);

    let result = kmer_map.into_hashmap(k_len);

    #[cfg(feature = "tracing")]
    info!(
        unique_kmers = result.len(),
        "K-mer counting with quality complete"
    );

    Ok(result)
}

/// Counts k-mers with progress callback.
///
/// Similar to [`count_kmers`], but invokes a callback after processing each sequence,
/// allowing callers to monitor progress during long-running operations.
/// Input format is auto-detected from the file extension.
///
/// # Arguments
///
/// * `path` - Path to the FASTA or FASTQ file
/// * `k` - K-mer length (must be 1-32)
/// * `callback` - Function called with a [`Progress`] snapshot after each sequence
///
/// # Returns
///
/// A `HashMap` mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be read
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::run::count_kmers_with_progress;
///
/// let counts = count_kmers_with_progress("genome.fa", 21, |progress| {
///     eprintln!(
///         "Processed {} sequences ({} bases)",
///         progress.sequences_processed,
///         progress.bases_processed
///     );
/// })?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn count_kmers_with_progress<P, F>(
    path: P,
    k: usize,
    callback: F,
) -> Result<HashMap<String, u64>, Box<dyn Error>>
where
    P: AsRef<Path> + Debug,
    F: Fn(Progress) + Send + Sync + 'static,
{
    use std::sync::Arc;

    #[cfg(feature = "tracing")]
    info!(k = k, path = ?path, "Starting k-mer counting with progress");

    // Validate k-mer length upfront to provide a clear error
    let k_len = KmerLength::new(k)?;

    #[cfg(feature = "tracing")]
    let read_span = info_span!("read_sequences", path = ?path).entered();

    let sequences = read(&path, SequenceFormat::Auto)?;

    #[cfg(feature = "tracing")]
    drop(read_span);

    #[cfg(feature = "tracing")]
    let process_span = info_span!("process_sequences").entered();

    let tracker = Arc::new(ProgressTracker::new());
    let callback = Arc::new(callback);
    let kmer_map = KmerMapWithProgress::new(tracker, callback).build(sequences, k);

    #[cfg(feature = "tracing")]
    drop(process_span);

    let result = kmer_map.into_hashmap(k_len);

    #[cfg(feature = "tracing")]
    info!(
        unique_kmers = result.len(),
        "K-mer counting with progress complete"
    );

    Ok(result)
}

/// Writes k-mer counts to stdout in the specified format.
///
/// # Arguments
///
/// * `counts` - `HashMap` of k-mer strings to their counts
/// * `format` - Output format (Fasta, Tsv, Json, or Histogram)
/// * `min_count` - Minimum count threshold (k-mers below this are excluded)
///
/// # Errors
///
/// Returns `ProcessError::WriteError` if output cannot be written.
/// Returns `ProcessError::JsonError` if JSON serialization fails.
#[allow(clippy::implicit_hasher)]
pub fn output_counts(
    counts: HashMap<String, u64>,
    format: OutputFormat,
    min_count: u64,
) -> Result<(), ProcessError> {
    let mut buf = BufWriter::new(stdout());
    let filtered: Vec<_> = counts
        .into_iter()
        .filter(|(_, count)| *count >= min_count)
        .collect();

    match format {
        OutputFormat::Fasta => {
            for (kmer, count) in filtered {
                writeln!(buf, ">{count}\n{kmer}")?;
            }
        }
        OutputFormat::Tsv => {
            for (kmer, count) in filtered {
                writeln!(buf, "{kmer}\t{count}")?;
            }
        }
        OutputFormat::Json => {
            let json_data: Vec<KmerCount> = filtered
                .into_iter()
                .map(|(kmer, count)| KmerCount { kmer, count })
                .collect();
            serde_json::to_writer_pretty(&mut buf, &json_data)?;
            writeln!(buf)?;
        }
        OutputFormat::Histogram => {
            use crate::histogram::compute_histogram;

            // Build a HashMap for compute_histogram
            let counts_map: HashMap<String, u64> = filtered.into_iter().collect();
            let histogram = compute_histogram(&counts_map);

            for (count, frequency) in histogram {
                writeln!(buf, "{count}\t{frequency}")?;
            }
        }
    }

    buf.flush()?;
    Ok(())
}

/// A custom [`DashMap`] w/ [`FxHasher`].
type DashFx = DashMap<u64, u64, BuildHasherDefault<FxHasher>>;

struct KmerMap(DashFx);

impl KmerMap {
    fn new() -> Self {
        Self(DashMap::with_hasher(
            BuildHasherDefault::<FxHasher>::default(),
        ))
    }

    fn build(self, sequences: rayon::vec::IntoIter<Bytes>, k: usize) -> Self {
        sequences.for_each(|seq| self.process_sequence(&seq, k));
        self
    }

    fn build_with_quality(
        self,
        sequences: rayon::vec::IntoIter<SequenceWithQuality>,
        k: usize,
        min_quality: Option<u8>,
    ) -> Self {
        sequences.for_each(|seq_qual| {
            self.process_sequence_with_quality(
                &seq_qual.seq,
                seq_qual.qual.as_deref(),
                k,
                min_quality,
            );
        });
        self
    }

    fn process_sequence(&self, seq: &Bytes, k: usize) {
        self.process_sequence_with_quality(seq, None, k, None);
    }

    fn process_sequence_with_quality(
        &self,
        seq: &Bytes,
        qual: Option<&[u8]>,
        k: usize,
        min_quality: Option<u8>,
    ) {
        if seq.len() < k {
            return;
        }

        // Pre-compute quality threshold as ASCII value (Phred+33)
        let quality_threshold = min_quality.map(|q| q.saturating_add(33));

        let mut i = 0;
        while i <= seq.len() - k {
            // Check quality if filtering is enabled
            if let (Some(q), Some(threshold)) = (qual, quality_threshold) {
                if let Some(bad_pos) = q[i..i + k].iter().position(|&qv| qv < threshold) {
                    i += bad_pos + 1; // Skip past low-quality base
                    continue;
                }
            }

            match pack_canonical(&seq[i..i + k], k) {
                Ok(canonical_bits) => {
                    self.0
                        .entry(canonical_bits)
                        .and_modify(|c| *c = c.saturating_add(1))
                        .or_insert(1);
                    i += 1;
                }
                Err(err) => {
                    i += err.position + 1;
                }
            }
        }
    }

    fn into_hashmap(self, k: KmerLength) -> HashMap<String, u64> {
        self.0
            .into_iter()
            .par_bridge()
            .map(|(packed_bits, count)| {
                let kmer_string = unpack_to_string(packed_bits, k);
                (kmer_string, count)
            })
            .collect()
    }
}

/// A k-mer map with progress tracking.
struct KmerMapWithProgress<F: Fn(Progress) + Send + Sync + 'static> {
    map: DashFx,
    tracker: std::sync::Arc<ProgressTracker>,
    callback: std::sync::Arc<F>,
}

impl<F: Fn(Progress) + Send + Sync + 'static> KmerMapWithProgress<F> {
    fn new(tracker: std::sync::Arc<ProgressTracker>, callback: std::sync::Arc<F>) -> Self {
        Self {
            map: DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default()),
            tracker,
            callback,
        }
    }

    #[allow(clippy::cast_possible_truncation)]
    fn build(self, sequences: rayon::vec::IntoIter<Bytes>, k: usize) -> Self {
        use rayon::prelude::ParallelIterator;

        sequences.for_each(|seq| {
            let len = seq.len() as u64;
            self.process_sequence(&seq, k);
            self.tracker.record_sequence(len);
            (self.callback)(self.tracker.snapshot());
        });
        self
    }

    fn process_sequence(&self, seq: &Bytes, k: usize) {
        if seq.len() < k {
            return;
        }
        let mut i = 0;

        while i <= seq.len() - k {
            match pack_canonical(&seq[i..i + k], k) {
                Ok(canonical_bits) => {
                    self.map
                        .entry(canonical_bits)
                        .and_modify(|c| *c = c.saturating_add(1))
                        .or_insert(1);
                    i += 1;
                }
                Err(err) => {
                    i += err.position + 1;
                }
            }
        }
    }

    fn into_hashmap(self, k: KmerLength) -> HashMap<String, u64> {
        self.map
            .into_iter()
            .par_bridge()
            .map(|(packed_bits, count)| {
                let kmer_string = unpack_to_string(packed_bits, k);
                (kmer_string, count)
            })
            .collect()
    }
}

/// Counts k-mers using memory-mapped I/O for potentially faster file access.
///
/// Memory-maps the FASTA file and processes it directly from the mapped region.
/// This can be more efficient for large files on systems with sufficient RAM.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A `HashMap` mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be opened or memory-mapped
/// - The file cannot be parsed as FASTA
///
/// # Safety
///
/// The underlying file must not be modified while being processed.
/// Modifying a mapped file leads to undefined behavior.
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::run::count_kmers_mmap;
///
/// let counts = count_kmers_mmap("large_genome.fa", 21)?;
/// println!("Found {} unique k-mers", counts.len());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[cfg(feature = "mmap")]
pub fn count_kmers_mmap<P>(path: P, k: usize) -> Result<HashMap<String, u64>, Box<dyn Error>>
where
    P: AsRef<Path> + Debug,
{
    use bio::io::fasta;
    use rayon::iter::IntoParallelIterator;
    use std::io::Cursor;

    #[cfg(feature = "tracing")]
    info!(k = k, path = ?path, "Starting memory-mapped k-mer counting");

    // Validate k-mer length upfront
    let k_len = KmerLength::new(k)?;

    #[cfg(feature = "tracing")]
    let mmap_span = info_span!("mmap_fasta", path = ?path).entered();

    let mmap =
        crate::mmap::MmapFasta::open(&path).map_err(|e| crate::error::KmeRustError::MmapError {
            source: e,
            path: path.as_ref().to_path_buf(),
        })?;

    #[cfg(feature = "tracing")]
    {
        drop(mmap_span);
        debug!(size_bytes = mmap.len(), "Memory-mapped file");
    }

    #[cfg(feature = "tracing")]
    let process_span = info_span!("process_sequences").entered();

    // Parse FASTA from the memory-mapped bytes
    let cursor = Cursor::new(mmap.as_bytes());
    let reader = fasta::Reader::new(cursor);
    let records: Vec<_> = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| crate::error::KmeRustError::SequenceParse {
            details: e.to_string(),
        })?;

    let kmer_map = KmerMap::new();
    let sequences: Vec<Bytes> = records
        .iter()
        .map(|r| Bytes::copy_from_slice(r.seq()))
        .collect();

    sequences
        .into_par_iter()
        .for_each(|seq| kmer_map.process_sequence(&seq, k));

    #[cfg(feature = "tracing")]
    drop(process_span);

    let result = kmer_map.into_hashmap(k_len);

    #[cfg(feature = "tracing")]
    info!(
        unique_kmers = result.len(),
        "Memory-mapped k-mer counting complete"
    );

    Ok(result)
}
