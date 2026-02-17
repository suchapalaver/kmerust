//! Property-based tests using proptest.
//!
//! These tests verify invariants that should hold across all valid inputs,
//! catching edge cases that might be missed by example-based tests.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::redundant_clone,
    clippy::explicit_iter_loop
)]

use bytes::Bytes;
use kmerust::index::{load_index, save_index, KmerIndex};
use kmerust::kmer::{pack_canonical, unpack_to_bytes, unpack_to_string, Kmer, KmerLength};
use kmerust::streaming::count_kmers_from_sequences;
use proptest::prelude::*;
use std::collections::HashMap;
use tempfile::NamedTempFile;

/// Strategy for generating valid DNA sequences of length 1-32.
fn dna_sequence(min_len: usize, max_len: usize) -> impl Strategy<Value = String> {
    proptest::collection::vec(
        prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')],
        min_len..=max_len,
    )
    .prop_map(|chars| chars.into_iter().collect())
}

/// Strategy for generating valid k-mer lengths (1-32).
fn kmer_length() -> impl Strategy<Value = usize> {
    1usize..=32
}

proptest! {
    /// Packing and unpacking a k-mer should be the identity operation.
    #[test]
    fn pack_unpack_roundtrip(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());
        let k = KmerLength::new(seq.len()).unwrap();

        let kmer = Kmer::from_sub(bytes).unwrap();
        let packed = kmer.pack();
        let unpacked = unpack_to_bytes(packed.packed_bits(), k);

        prop_assert_eq!(unpacked.as_ref(), seq.as_bytes());
    }

    /// Unpacking to string should produce the same result as unpacking to bytes.
    #[test]
    fn unpack_bytes_equals_string(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());
        let k = KmerLength::new(seq.len()).unwrap();

        let kmer = Kmer::from_sub(bytes).unwrap();
        let packed = kmer.pack();

        let unpacked_bytes = unpack_to_bytes(packed.packed_bits(), k);
        let unpacked_string = unpack_to_string(packed.packed_bits(), k);

        prop_assert_eq!(unpacked_bytes.as_ref(), unpacked_string.as_bytes());
    }

    /// The canonical form should be idempotent: canonical(canonical(x)) == canonical(x).
    #[test]
    fn canonical_is_idempotent(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());

        let canonical1 = Kmer::from_sub(bytes.clone()).unwrap().pack().canonical();
        let canonical_bytes = canonical1.bytes().clone();

        // Taking canonical of the already-canonical form should give same result
        let canonical2 = Kmer::from_sub(canonical_bytes).unwrap().pack().canonical();

        prop_assert_eq!(canonical1.packed_bits(), canonical2.packed_bits());
    }

    /// A k-mer and its reverse complement should have the same canonical form.
    #[test]
    fn kmer_and_rc_have_same_canonical(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());

        // Compute reverse complement manually
        let rc: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => unreachable!(),
            })
            .collect();

        let canonical1 = Kmer::from_sub(bytes).unwrap().pack().canonical();
        let canonical2 = Kmer::from_sub(Bytes::from(rc)).unwrap().pack().canonical();

        prop_assert_eq!(canonical1.packed_bits(), canonical2.packed_bits());
    }

    /// The canonical form should be lexicographically <= both the original and its RC.
    #[test]
    fn canonical_is_lexicographically_smallest(seq in dna_sequence(1, 32)) {
        let bytes = Bytes::from(seq.clone());

        // Compute reverse complement
        let rc: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => unreachable!(),
            })
            .collect();

        let canonical = Kmer::from_sub(bytes).unwrap().pack().canonical();
        let canonical_str = std::str::from_utf8(canonical.bytes()).unwrap();

        // Canonical should be <= both original and reverse complement
        prop_assert!(canonical_str <= seq.as_str());
        prop_assert!(canonical_str <= rc.as_str());
    }

    /// KmerLength should accept values 1-32 and reject others.
    #[test]
    fn kmer_length_accepts_valid_range(k in 1usize..=32) {
        let result = KmerLength::new(k);
        prop_assert!(result.is_ok());
        prop_assert_eq!(result.unwrap().get(), k);
    }

    #[test]
    fn kmer_length_rejects_zero(_dummy in Just(())) {
        let result = KmerLength::new(0);
        prop_assert!(result.is_err());
    }

    #[test]
    fn kmer_length_rejects_too_large(k in 33usize..1000) {
        let result = KmerLength::new(k);
        prop_assert!(result.is_err());
    }

    /// Soft-masked (lowercase) sequences should produce same result as uppercase.
    #[test]
    fn soft_masked_equals_uppercase(seq in dna_sequence(1, 32)) {
        let uppercase = Bytes::from(seq.clone());
        let lowercase = Bytes::from(seq.to_lowercase());

        let packed_upper = Kmer::from_sub(uppercase).unwrap().pack();
        let packed_lower = Kmer::from_sub(lowercase).unwrap().pack();

        prop_assert_eq!(packed_upper.packed_bits(), packed_lower.packed_bits());
    }

    /// Mixed case should produce same result as all uppercase.
    #[test]
    fn mixed_case_equals_uppercase(seq in dna_sequence(1, 32)) {
        let uppercase = Bytes::from(seq.clone());

        // Create mixed case version
        let mixed: String = seq
            .chars()
            .enumerate()
            .map(|(i, c)| if i % 2 == 0 { c } else { c.to_ascii_lowercase() })
            .collect();
        let mixed_bytes = Bytes::from(mixed);

        let packed_upper = Kmer::from_sub(uppercase).unwrap().pack();
        let packed_mixed = Kmer::from_sub(mixed_bytes).unwrap().pack();

        prop_assert_eq!(packed_upper.packed_bits(), packed_mixed.packed_bits());
    }

    /// Packing should be deterministic.
    #[test]
    fn packing_is_deterministic(seq in dna_sequence(1, 32)) {
        let bytes1 = Bytes::from(seq.clone());
        let bytes2 = Bytes::from(seq);

        let packed1 = Kmer::from_sub(bytes1).unwrap().pack();
        let packed2 = Kmer::from_sub(bytes2).unwrap().pack();

        prop_assert_eq!(packed1.packed_bits(), packed2.packed_bits());
    }

    /// Different sequences of the same length should have different packed representations.
    /// This tests that the packing is injective.
    #[test]
    fn different_sequences_different_packing(
        seq in dna_sequence(1, 16),
        mutation_pos in 0usize..16,
    ) {
        // Skip if mutation position is out of bounds
        prop_assume!(mutation_pos < seq.len());

        // Create a second sequence by changing one base
        let original_char = seq.chars().nth(mutation_pos).unwrap();
        let new_char = match original_char {
            'A' => 'C',
            'C' => 'G',
            'G' => 'T',
            'T' => 'A',
            _ => unreachable!(),
        };

        let mut chars: Vec<char> = seq.chars().collect();
        chars[mutation_pos] = new_char;
        let seq2: String = chars.into_iter().collect();

        let packed1 = Kmer::from_sub(Bytes::from(seq.clone())).unwrap().pack();
        let packed2 = Kmer::from_sub(Bytes::from(seq2)).unwrap().pack();

        prop_assert_ne!(packed1.packed_bits(), packed2.packed_bits());
    }

    /// The length of unpacked bytes should equal the original k-mer length.
    #[test]
    fn unpack_preserves_length(k in kmer_length(), bits in any::<u64>()) {
        let k_len = KmerLength::new(k).unwrap();
        let unpacked = unpack_to_bytes(bits, k_len);

        prop_assert_eq!(unpacked.len(), k);
    }

    /// Unpacked bytes should only contain valid bases (A, C, G, T).
    #[test]
    fn unpack_produces_valid_bases(k in kmer_length(), bits in any::<u64>()) {
        let k_len = KmerLength::new(k).unwrap();
        let unpacked = unpack_to_bytes(bits, k_len);

        for &byte in unpacked.iter() {
            prop_assert!(
                byte == b'A' || byte == b'C' || byte == b'G' || byte == b'T',
                "Invalid base: {}", byte as char
            );
        }
    }

    /// Index save/load roundtrip should preserve all entries exactly.
    ///
    /// Property: load(save(index)) = index
    #[test]
    fn index_roundtrip_preserves_all_entries(
        k in 1usize..=32,
        entries in proptest::collection::hash_map(any::<u64>(), 1u64..1000, 0..100)
    ) {
        let k_len = KmerLength::new(k).unwrap();
        let index = KmerIndex::new(k_len, entries.clone());

        let tmp = NamedTempFile::with_suffix(".kmix").unwrap();
        save_index(&index, tmp.path()).unwrap();
        let loaded = load_index(tmp.path()).unwrap();

        prop_assert_eq!(loaded.k(), k_len);
        prop_assert_eq!(loaded.counts(), &entries);
    }

    /// Total k-mer count should not exceed the number of valid windows in the input.
    ///
    /// Property: Σ(counts) ≤ (seq.len - k + 1) for a single sequence
    #[test]
    fn total_count_at_most_valid_windows(
        seq in dna_sequence(1, 100),
        k in 1usize..=32
    ) {
        prop_assume!(seq.len() >= k);

        let k_len = KmerLength::new(k).unwrap();
        let counts = count_kmers_from_sequences(
            vec![Bytes::from(seq.clone())].into_iter(),
            k_len
        );

        let total: u64 = counts.values().sum();
        let max_windows = (seq.len() - k + 1) as u64;

        prop_assert!(
            total <= max_windows,
            "Total count {total} exceeds max windows {max_windows}"
        );
    }

    /// A k-mer and its reverse complement should be counted together under one canonical entry.
    ///
    /// Property: When counting a k-mer and its reverse complement as separate sequences
    /// (each of exactly length k), they should produce exactly one canonical entry
    /// with count 2 (or count 2 for palindromes since they're the same sequence).
    #[test]
    fn kmer_and_rc_count_together(seq in dna_sequence(1, 32)) {
        // Compute reverse complement
        let rc: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => unreachable!(),
            })
            .collect();

        let k = seq.len();
        let k_len = KmerLength::new(k).unwrap();

        // Count the k-mer and its RC as separate sequences (each exactly k bases)
        // This ensures we get exactly one k-mer from each sequence
        let counts: HashMap<u64, u64> = count_kmers_from_sequences(
            vec![Bytes::from(seq.clone()), Bytes::from(rc.clone())].into_iter(),
            k_len
        );

        // Both should map to the same canonical form
        prop_assert_eq!(
            counts.len(), 1,
            "K-mer and RC should produce exactly 1 canonical entry, got {}", counts.len()
        );

        // Count should be 2 (once from original, once from RC)
        let kmer_count = *counts.values().next().unwrap();
        prop_assert_eq!(
            kmer_count, 2,
            "K-mer and RC should have combined count 2, got {}", kmer_count
        );
    }

    /// The zero-allocation fast path should produce identical results to the
    /// type-state pipeline for all valid DNA sequences.
    ///
    /// Property: pack_canonical(seq) == Kmer::from_sub(seq).pack().canonical().packed_bits()
    #[test]
    fn fast_path_matches_type_state_pipeline(seq in dna_sequence(1, 32)) {
        let via_pipeline = Kmer::from_sub(Bytes::from(seq.clone()))
            .unwrap()
            .pack()
            .canonical()
            .packed_bits();
        let via_fast_path = pack_canonical(seq.as_bytes()).unwrap();

        prop_assert_eq!(
            via_pipeline, via_fast_path,
            "mismatch for sequence {}", seq
        );
    }
}
