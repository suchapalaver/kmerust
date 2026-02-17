# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.2] - 2026-02-17

### Performance

- 2.4-10x faster k-mer counting by eliminating per-window heap allocations

### Fixed

- `validate_and_pack` now rejects empty input and sequences longer than 32 in debug builds
- `pack_canonical` catches empty slices in debug builds
- Panic tests for empty slices pass correctly under release builds
- `reverse_complement_bits` documents release-mode behaviour for out-of-range k values

### Internal

- Simplified `pack_canonical` API and removed dead code
- Added property tests for fast-path correctness and invalid k-mer length detection
- Added debug assertions and safety comments on public bit-packing functions

## [0.3.1] - 2025-01-27

### Internal

- Enabled pedantic, nursery, and all clippy lint groups in `Cargo.toml` for comprehensive static analysis
- Added property tests for index serialization roundtrip and canonical k-mer deduplication
- Made numerous functions `const` where applicable
- Improved code quality with explicit safety annotations on intentional `unsafe` blocks

## [0.3.0] - 2026-01-27

### Added

#### FASTQ Support

- Full FASTQ file format support with automatic detection from file extension (`.fq`, `.fastq`, `.fq.gz`, `.fastq.gz`)
- `--input-format` / `-i` CLI flag for explicit format specification (required for stdin FASTQ)
- `SequenceFormat` enum (`Auto`, `Fasta`, `Fastq`) for programmatic format selection
- FASTQ support in streaming APIs (`count_kmers_streaming`, `count_kmers_sequential`)
- Gzip-compressed FASTQ support (`.fq.gz`)

#### Quality-Based Filtering

- `--min-quality` / `-Q` CLI flag to filter k-mers by base quality scores (Phred 0-93)
- K-mers containing bases below the quality threshold are skipped
- Smart skip optimization: jumps past low-quality bases rather than checking overlapping windows
- `count_kmers_with_quality` function for programmatic quality filtering
- Appropriate warnings when quality filtering is used with FASTA or stdin

#### Histogram Output

- `--format histogram` output mode for k-mer frequency spectrum (count of counts)
- Tab-separated output: `count<TAB>frequency`, sorted by count ascending
- `KmerHistogram` type alias (`BTreeMap<u64, u64>`) for histogram data
- `compute_histogram` and `compute_histogram_packed` functions
- `histogram_stats` function returning `HistogramStats` (total, distinct, mode, mean)
- `KmerCounter::histogram()` builder method

#### K-mer Index Serialization

- `--save <PATH>` flag to save k-mer counts to binary index file (`.kmix`)
- `kmerust query <INDEX> <KMER>` subcommand to look up k-mer counts from saved index
- Compact binary format with versioned header and CRC32 integrity checking
- Automatic gzip compression when saving to `.kmix.gz`
- `KmerIndex` struct with `save_index` and `load_index` functions
- Query support: case-insensitive, canonicalized lookups

### Changed

- Package description updated to mention FASTQ support
- Error types renamed for format neutrality: `FastaRead` → `SequenceRead`, `FastaParse` → `SequenceParse`
- Consolidated duplicate output code in main.rs for better error handling
- `output_counts` function is now public in the `run` module

### Fixed

- Quality filtering correctly warns when used with stdin (not yet supported there)
- Quality filtering correctly warns when used with FASTA input (no quality data available)

## [0.2.1] - 2026-01-23

### Added

- Stdin input support for Unix pipeline integration - read FASTA data from stdin using `-` or by omitting the path argument
- `count_kmers_stdin` function for programmatic stdin reading
- `count_kmers_from_reader` function for counting k-mers from any `BufRead` source

### Changed

- CI updated to actions/checkout v6

## [0.2.0] - 2026-01-22

### Added

#### Type Safety & API Correctness

- `KmerLength` newtype that enforces valid k-mer lengths (1-32) at construction time
- Type-state pattern for `Kmer<S>` with `Unpacked`, `Packed`, and `Canonical` states
- Exhaustive error types (`KmeRustError`, `KmerLengthError`, `InvalidBaseError`) replacing `Box<dyn Error>`
- Compile-time enforcement of k-mer operation ordering

#### Performance Hardening

- Compile-time lookup tables (`PACK_TABLE`, `COMPLEMENT_TABLE`, `UNPACK_TABLE`) for bit packing
- Zero-allocation canonical comparison using lazy iterator comparison
- Streaming API (`count_kmers_streaming`, `count_kmers_streaming_packed`) for memory-efficient processing
- `count_kmers_from_sequences` for counting from in-memory byte slices

#### Testing Rigor

- Direct library API test suite (`tests/library_tests.rs`)
- Property-based tests using proptest (`tests/property_tests.rs`)
- Jellyfish compatibility tests (`tests/jellyfish_compat.rs`)
- Fuzz targets for input handling (`fuzz/fuzz_targets/`)

#### Production Readiness

- Gzip compressed input support (`gzip` feature)
- Memory-mapped I/O for large files (`mmap` feature)
- Progress callback API (`count_kmers_with_progress`)
- Tracing instrumentation (`tracing` feature)

#### API Polish

- Builder pattern API (`KmerCounter`) for ergonomic k-mer counting
- Async API with Tokio support (`async` feature)
  - `count_kmers_async` for non-blocking k-mer counting
  - `AsyncKmerCounter` builder for async operations
- Example programs demonstrating common use cases:
  - `basic_count.rs` - Simple k-mer counting
  - `streaming_large_file.rs` - Memory-efficient processing
  - `compare_with_jellyfish.rs` - Validation against Jellyfish
  - `progress_bar.rs` - Progress reporting during counting
  - `async_count.rs` - Async k-mer counting
- `full` feature flag enabling all optional features

#### Infrastructure

- Comprehensive test suite with unit and integration tests
- CI/CD pipeline with GitHub Actions
- Rustdoc documentation for public APIs
- MSRV (Minimum Supported Rust Version) set to 1.75
- Semver checking with cargo-semver-checks in CI

### Changed

- `count_kmers` now validates k-mer length upfront and returns a clear error
- K-mer packing uses direct lookup tables instead of match statements
- Canonical computation avoids allocation when the original form is already canonical
- Simplified `process_valid_kmer()` by removing misleading "optimization" that only helped ~50% of k-mers
- Made `needletail` a truly optional dependency (compile-time savings when not using needletail feature)
- Builder API always uses canonical k-mers (removed unused `canonical` field)
- Updated all dependencies to latest stable versions
- Replaced `fxhash` with `rustc-hash` (maintained successor)
- CLI version now reads from Cargo.toml
- Improved error handling in FASTA readers

### Removed

- Deprecated `config.rs` module
- Unused `canonical` field and method from `KmerCounter` builder (all counting is canonical)

### Fixed

- Panic when sequence length is smaller than k-mer length
- Clippy warnings

## [0.1.0] - Initial Release

### Added

- Parallel k-mer counting using rayon and dashmap
- Canonical k-mer computation
- Support for rust-bio and needletail FASTA readers
- K-mer lengths from 1 to 32
- Jellyfish-compatible output format

[Unreleased]: https://github.com/suchapalaver/kmerust/compare/v0.3.2...HEAD
[0.3.2]: https://github.com/suchapalaver/kmerust/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/suchapalaver/kmerust/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/suchapalaver/kmerust/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/suchapalaver/kmerust/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/suchapalaver/kmerust/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/suchapalaver/kmerust/releases/tag/v0.1.0
