#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::uninlined_format_args,
    clippy::semicolon_if_nothing_returned
)]

use bytes::Bytes;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use kmerust::kmer::{
    pack_canonical, reverse_complement_bits, unpack_to_bytes, validate_and_pack, Kmer, KmerLength,
};
use kmerust::run::count_kmers;
use kmerust::streaming::{count_kmers_from_sequences, count_kmers_streaming};
use std::io::Write;
use tempfile::NamedTempFile;

fn bench_from_sub(c: &mut Criterion) {
    let mut group = c.benchmark_group("Kmer::from_sub");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| Kmer::from_sub(black_box(bytes.clone())))
        });
    }

    group.finish();
}

fn bench_pack(c: &mut Criterion) {
    let mut group = c.benchmark_group("Kmer::pack");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| {
                let kmer = Kmer::from_sub(bytes.clone()).unwrap();
                black_box(kmer.pack())
            })
        });
    }

    group.finish();
}

fn bench_canonical(c: &mut Criterion) {
    let mut group = c.benchmark_group("Kmer::canonical");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| {
                let kmer = Kmer::from_sub(bytes.clone()).unwrap();
                black_box(kmer.pack().canonical())
            })
        });
    }

    group.finish();
}

fn bench_canonical_no_alloc(c: &mut Criterion) {
    // Benchmark canonical when original is already smallest (no allocation needed)
    let mut group = c.benchmark_group("Kmer::canonical_no_alloc");

    for k in [5, 11, 21, 31] {
        // "AAAA..." is already canonical (smaller than "TTTT...")
        let seq = "A".repeat(k);
        let bytes = Bytes::copy_from_slice(seq.as_bytes());

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| {
                let kmer = Kmer::from_sub(bytes.clone()).unwrap();
                black_box(kmer.pack().canonical())
            })
        });
    }

    group.finish();
}

fn bench_canonical_needs_alloc(c: &mut Criterion) {
    // Benchmark canonical when reverse complement is smaller (allocation needed)
    let mut group = c.benchmark_group("Kmer::canonical_needs_alloc");

    for k in [5, 11, 21, 31] {
        // "TTTT..." needs allocation for reverse complement "AAAA..."
        let seq = "T".repeat(k);
        let bytes = Bytes::copy_from_slice(seq.as_bytes());

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| {
                let kmer = Kmer::from_sub(bytes.clone()).unwrap();
                black_box(kmer.pack().canonical())
            })
        });
    }

    group.finish();
}

fn bench_unpack(c: &mut Criterion) {
    let mut group = c.benchmark_group("unpack_to_bytes");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);
        let packed = Kmer::from_sub(bytes).unwrap().pack();
        let packed_bits = packed.packed_bits();
        let k_len = KmerLength::new(k).unwrap();

        group.bench_with_input(
            BenchmarkId::from_parameter(k),
            &(packed_bits, k_len),
            |b, &(bits, k_len)| b.iter(|| black_box(unpack_to_bytes(bits, k_len))),
        );
    }

    group.finish();
}

fn bench_count_kmers_small(c: &mut Criterion) {
    let mut group = c.benchmark_group("count_kmers");

    // Create a small test file
    let mut file = NamedTempFile::new().unwrap();
    for i in 0..100 {
        writeln!(file, ">seq{i}").unwrap();
        writeln!(file, "{}", "ACGTACGTACGTACGTACGTACGTACGTACGT".repeat(10)).unwrap();
    }
    let path = file.path().to_path_buf();

    for k in [5, 11, 21] {
        group.bench_with_input(BenchmarkId::from_parameter(k), &k, |b, &k| {
            b.iter(|| count_kmers(black_box(&path), black_box(k)))
        });
    }

    group.finish();
}

fn bench_count_kmers_streaming(c: &mut Criterion) {
    let mut group = c.benchmark_group("count_kmers_streaming");

    // Create a small test file
    let mut file = NamedTempFile::new().unwrap();
    for i in 0..100 {
        writeln!(file, ">seq{i}").unwrap();
        writeln!(file, "{}", "ACGTACGTACGTACGTACGTACGTACGTACGT".repeat(10)).unwrap();
    }
    let path = file.path().to_path_buf();

    for k in [5, 11, 21] {
        group.bench_with_input(BenchmarkId::from_parameter(k), &k, |b, &k| {
            b.iter(|| count_kmers_streaming(black_box(&path), black_box(k)))
        });
    }

    group.finish();
}

fn bench_count_from_sequences(c: &mut Criterion) {
    let mut group = c.benchmark_group("count_kmers_from_sequences");

    // Pre-create sequences in memory
    let sequences: Vec<Bytes> = (0..100)
        .map(|_| Bytes::from("ACGTACGTACGTACGTACGTACGTACGTACGT".repeat(10)))
        .collect();

    for k in [5, 11, 21] {
        let k_len = KmerLength::new(k).unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(k), &k_len, |b, &k_len| {
            b.iter(|| {
                count_kmers_from_sequences(
                    black_box(sequences.clone().into_iter()),
                    black_box(k_len),
                )
            })
        });
    }

    group.finish();
}

fn bench_validate_and_pack(c: &mut Criterion) {
    let mut group = c.benchmark_group("validate_and_pack");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = &seq.as_bytes()[..k];

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| validate_and_pack(black_box(bytes)))
        });
    }

    group.finish();
}

fn bench_pack_canonical(c: &mut Criterion) {
    let mut group = c.benchmark_group("pack_canonical");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = &seq.as_bytes()[..k];

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| pack_canonical(black_box(bytes)))
        });
    }

    group.finish();
}

fn bench_pack_canonical_needs_rc(c: &mut Criterion) {
    let mut group = c.benchmark_group("pack_canonical_needs_rc");

    for k in [5, 11, 21, 31] {
        let seq = "T".repeat(k);
        let bytes = seq.as_bytes();

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| pack_canonical(black_box(bytes)))
        });
    }

    group.finish();
}

fn bench_reverse_complement_bits(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse_complement_bits");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);
        let packed = Kmer::from_sub(bytes).unwrap().pack();
        let bits = packed.packed_bits();

        group.bench_with_input(
            BenchmarkId::from_parameter(k),
            &(bits, k),
            |b, &(bits, k)| b.iter(|| reverse_complement_bits(black_box(bits), black_box(k))),
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_from_sub,
    bench_pack,
    bench_canonical,
    bench_canonical_no_alloc,
    bench_canonical_needs_alloc,
    bench_unpack,
    bench_validate_and_pack,
    bench_pack_canonical,
    bench_pack_canonical_needs_rc,
    bench_reverse_complement_bits,
    bench_count_kmers_small,
    bench_count_kmers_streaming,
    bench_count_from_sequences,
);

criterion_main!(benches);
