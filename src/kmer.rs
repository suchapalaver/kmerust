//! K-mer representation and manipulation.
//!
//! This module provides types for working with DNA k-mers, including:
//! - Bit-packing k-mers into 64-bit integers (supports k <= 32)
//! - Computing canonical k-mers (lexicographically smaller of k-mer/reverse-complement)
//! - Converting between byte and packed representations

use bytes::Bytes;

use crate::error::KmerLengthError;

// ============================================================================
// Compile-time Lookup Tables for Performance
// ============================================================================

/// Lookup table for packing DNA bases into 2-bit representation.
/// Maps ASCII byte values directly to their 2-bit encoding: A=0, C=1, G=2, T=3.
/// Invalid bytes map to 0 (but should never be looked up after validation).
const PACK_TABLE: [u64; 256] = {
    let mut table = [0u64; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// Lookup table for unpacking 2-bit values back to ASCII bytes.
const UNPACK_TABLE: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// A validated k-mer length (1-32 inclusive).
///
/// This newtype enforces that k-mer lengths are within the valid range for
/// 64-bit packed representation (2 bits per base, so max 32 bases).
///
/// # Invariant
///
/// The inner value is always in [`Self::MIN`]..=[`Self::MAX`], enforced at construction.
///
/// # Examples
///
/// ```
/// use kmerust::kmer::KmerLength;
///
/// let k = KmerLength::new(21).unwrap();
/// assert_eq!(k.get(), 21);
///
/// // Invalid lengths are rejected
/// assert!(KmerLength::new(0).is_err());
/// assert!(KmerLength::new(33).is_err());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct KmerLength(u8);

impl KmerLength {
    /// Minimum valid k-mer length.
    pub const MIN: u8 = 1;

    /// Maximum valid k-mer length (limited by 64-bit packing: 2 bits × 32 = 64 bits).
    pub const MAX: u8 = 32;

    /// Creates a new `KmerLength` from a `usize` value.
    ///
    /// # Errors
    ///
    /// Returns [`KmerLengthError`] if `k` is outside the valid range (1-32).
    ///
    /// # Examples
    ///
    /// ```
    /// use kmerust::kmer::KmerLength;
    ///
    /// assert!(KmerLength::new(1).is_ok());
    /// assert!(KmerLength::new(32).is_ok());
    /// assert!(KmerLength::new(0).is_err());
    /// assert!(KmerLength::new(33).is_err());
    /// ```
    #[allow(clippy::cast_possible_truncation)]
    pub const fn new(k: usize) -> Result<Self, KmerLengthError> {
        if k < Self::MIN as usize || k > Self::MAX as usize {
            return Err(KmerLengthError {
                k,
                min: Self::MIN,
                max: Self::MAX,
            });
        }
        // SAFETY: range check above guarantees k fits in u8
        Ok(Self(k as u8))
    }

    /// Creates a `KmerLength` without validation.
    ///
    /// # Safety
    ///
    /// The caller must ensure that `k` is in the range 1..=32.
    /// Using a value outside this range will lead to incorrect bit-packing.
    ///
    /// # Examples
    ///
    /// ```
    /// use kmerust::kmer::KmerLength;
    ///
    /// // SAFETY: 21 is within the valid range 1..=32
    /// let k = unsafe { KmerLength::new_unchecked(21) };
    /// assert_eq!(k.get(), 21);
    /// ```
    #[inline]
    #[allow(unsafe_code)]
    pub const unsafe fn new_unchecked(k: u8) -> Self {
        Self(k)
    }

    /// Returns the k-mer length as a `usize`.
    #[inline]
    pub const fn get(self) -> usize {
        self.0 as usize
    }

    /// Returns the k-mer length as a `u8`.
    #[inline]
    pub const fn as_u8(self) -> u8 {
        self.0
    }
}

impl TryFrom<usize> for KmerLength {
    type Error = KmerLengthError;

    fn try_from(k: usize) -> Result<Self, Self::Error> {
        Self::new(k)
    }
}

impl From<KmerLength> for usize {
    fn from(k: KmerLength) -> Self {
        k.get()
    }
}

use std::marker::PhantomData;

use crate::error::InvalidBaseError;

// ============================================================================
// Type State Markers
// ============================================================================

/// Marker type: k-mer has validated bytes but is not yet packed.
#[derive(Debug, Clone, Copy, Default)]
pub struct Unpacked;

/// Marker type: k-mer bytes have been packed into a 64-bit integer.
#[derive(Debug, Clone, Copy, Default)]
pub struct Packed;

/// Marker type: k-mer is in canonical form (lexicographically smaller of k-mer/RC).
#[derive(Debug, Clone, Copy, Default)]
pub struct Canonical;

// ============================================================================
// Kmer Type with Type States
// ============================================================================

/// A DNA k-mer with compile-time state tracking.
///
/// K-mers progress through states via method calls:
/// - [`Kmer<Unpacked>`]: Created from validated bytes
/// - [`Kmer<Packed>`]: Bytes packed into 64-bit representation
/// - [`Kmer<Canonical>`]: Converted to canonical form
///
/// # Type Safety
///
/// The type state pattern ensures at compile time that operations are performed
/// in the correct order. For example, you cannot access `packed_bits()` on an
/// unpacked k-mer.
///
/// # Examples
///
/// ```
/// use bytes::Bytes;
/// use kmerust::kmer::Kmer;
///
/// let unpacked = Kmer::from_sub(Bytes::from_static(b"GATTACA")).unwrap();
/// let packed = unpacked.pack();
/// let canonical = packed.canonical();
///
/// // Access the canonical packed bits for use as a hash key
/// let bits = canonical.packed_bits();
/// ```
#[derive(Debug, Clone)]
pub struct Kmer<S = Unpacked> {
    bytes: Bytes,
    packed_bits: u64,
    is_reverse_complement: bool,
    _state: PhantomData<S>,
}

impl Default for Kmer<Unpacked> {
    fn default() -> Self {
        Self {
            bytes: Bytes::new(),
            packed_bits: 0,
            is_reverse_complement: false,
            _state: PhantomData,
        }
    }
}

impl<S> Kmer<S> {
    /// Returns the k-mer bytes.
    #[inline]
    pub const fn bytes(&self) -> &Bytes {
        &self.bytes
    }
}

impl Kmer<Unpacked> {
    /// Creates a k-mer from a byte sequence.
    ///
    /// Validates that all bytes are valid DNA bases (A, C, G, T, or lowercase).
    /// Soft-masked (lowercase) bases are converted to uppercase.
    ///
    /// # Errors
    ///
    /// Returns [`InvalidBaseError`] with the position and value of the first invalid byte.
    ///
    /// # Examples
    ///
    /// ```
    /// use bytes::Bytes;
    /// use kmerust::kmer::Kmer;
    ///
    /// let kmer = Kmer::from_sub(Bytes::from_static(b"GATTACA")).unwrap();
    /// assert_eq!(kmer.bytes().as_ref(), b"GATTACA");
    ///
    /// // Lowercase bases are normalized
    /// let kmer = Kmer::from_sub(Bytes::from_static(b"gattaca")).unwrap();
    /// assert_eq!(kmer.bytes().as_ref(), b"GATTACA");
    ///
    /// // Invalid bases return an error
    /// let result = Kmer::from_sub(Bytes::from_static(b"GANTACA"));
    /// assert!(result.is_err());
    /// ```
    #[allow(clippy::needless_pass_by_value)]
    pub fn from_sub(sub: Bytes) -> Result<Self, InvalidBaseError> {
        let normalized: Result<Vec<u8>, InvalidBaseError> = sub
            .iter()
            .enumerate()
            .map(|(i, &byte)| match byte {
                b'A' | b'C' | b'G' | b'T' => Ok(byte),
                b'a' | b'c' | b'g' | b't' => Ok(byte.to_ascii_uppercase()),
                _ => Err(InvalidBaseError {
                    base: byte,
                    position: i,
                }),
            })
            .collect();

        Ok(Self {
            bytes: Bytes::from(normalized?),
            packed_bits: 0,
            is_reverse_complement: false,
            _state: PhantomData,
        })
    }

    /// Packs the k-mer bytes into a 64-bit integer.
    ///
    /// Each base is encoded as 2 bits: A=00, C=01, G=10, T=11.
    /// Consumes the unpacked k-mer and returns a packed k-mer.
    ///
    /// # Examples
    ///
    /// ```
    /// use bytes::Bytes;
    /// use kmerust::kmer::Kmer;
    ///
    /// let unpacked = Kmer::from_sub(Bytes::from_static(b"ACGT")).unwrap();
    /// let packed = unpacked.pack();
    /// // ACGT = 00 01 10 11 = 0b00011011 = 27
    /// assert_eq!(packed.packed_bits(), 0b00_01_10_11);
    /// ```
    pub fn pack(self) -> Kmer<Packed> {
        let packed_bits = pack_bytes(&self.bytes);
        Kmer {
            bytes: self.bytes,
            packed_bits,
            is_reverse_complement: false,
            _state: PhantomData,
        }
    }
}

impl Kmer<Packed> {
    /// Returns the packed 64-bit representation.
    #[inline]
    pub const fn packed_bits(&self) -> u64 {
        self.packed_bits
    }

    /// Converts the k-mer to its canonical form.
    ///
    /// The canonical form is the lexicographically smaller of the k-mer
    /// and its reverse complement. This ensures that a k-mer and its
    /// reverse complement are treated as equivalent.
    ///
    /// # Examples
    ///
    /// ```
    /// use bytes::Bytes;
    /// use kmerust::kmer::Kmer;
    ///
    /// // TTT -> AAA (reverse complement), AAA is smaller
    /// let kmer = Kmer::from_sub(Bytes::from_static(b"TTT")).unwrap()
    ///     .pack()
    ///     .canonical();
    /// assert_eq!(kmer.bytes().as_ref(), b"AAA");
    /// assert!(kmer.is_reverse_complement());
    ///
    /// // AAA is already canonical
    /// let kmer = Kmer::from_sub(Bytes::from_static(b"AAA")).unwrap()
    ///     .pack()
    ///     .canonical();
    /// assert_eq!(kmer.bytes().as_ref(), b"AAA");
    /// assert!(!kmer.is_reverse_complement());
    /// ```
    pub fn canonical(self) -> Kmer<Canonical> {
        let k = self.bytes.len();
        let rc_bits = reverse_complement_bits(self.packed_bits, k);
        let use_reverse_complement = rc_bits < self.packed_bits;

        if use_reverse_complement {
            // Reconstruct RC bytes from packed bits (for public API consumers)
            let rc_bytes = unpack_to_bytes_from_raw(rc_bits, k);
            Kmer {
                bytes: rc_bytes,
                packed_bits: rc_bits,
                is_reverse_complement: true,
                _state: PhantomData,
            }
        } else {
            Kmer {
                bytes: self.bytes,
                packed_bits: self.packed_bits,
                is_reverse_complement: false,
                _state: PhantomData,
            }
        }
    }
}

impl Kmer<Canonical> {
    /// Returns the canonical packed 64-bit representation.
    #[inline]
    pub const fn packed_bits(&self) -> u64 {
        self.packed_bits
    }

    /// Returns whether this k-mer was converted to its reverse complement.
    #[inline]
    pub const fn is_reverse_complement(&self) -> bool {
        self.is_reverse_complement
    }
}

// ============================================================================
// Unpacking: Create Kmer from packed bits
// ============================================================================

/// Unpacks a 64-bit integer into k-mer bytes.
///
/// This is the inverse of packing: given packed bits and a k-mer length,
/// reconstructs the DNA sequence.
///
/// # Arguments
///
/// * `packed_bits` - The packed 64-bit representation
/// * `k` - The k-mer length
///
/// # Examples
///
/// ```
/// use kmerust::kmer::{unpack_to_bytes, KmerLength};
///
/// let k = KmerLength::new(4).unwrap();
/// // 0b00011011 = ACGT
/// let bytes = unpack_to_bytes(0b00_01_10_11, k);
/// assert_eq!(bytes.as_ref(), b"ACGT");
/// ```
pub fn unpack_to_bytes(packed_bits: u64, k: KmerLength) -> Bytes {
    let k = k.get();
    (0..k)
        .map(|i| {
            let shift = (k - 1 - i) * 2;
            let bits = ((packed_bits >> shift) & 0b11) as usize;
            UNPACK_TABLE[bits]
        })
        .collect()
}

/// Unpacks packed bits to a String.
///
/// Convenience function for when a String representation is needed.
///
/// # Safety Note
///
/// This function uses `String::from_utf8_unchecked` internally because
/// the unpacking process only produces valid ASCII bytes (A, C, G, T).
#[allow(unsafe_code)]
pub fn unpack_to_string(packed_bits: u64, k: KmerLength) -> String {
    let bytes = unpack_to_bytes(packed_bits, k);
    // SAFETY: unpack_to_bytes only produces bytes from KmerByte (A, C, G, T),
    // which are valid ASCII and therefore valid UTF-8
    unsafe { String::from_utf8_unchecked(bytes.to_vec()) }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Packs DNA bytes into a 64-bit integer using lookup table.
///
/// Each base is encoded as 2 bits: A=00, C=01, G=10, T=11.
/// Uses a compile-time lookup table for direct indexing instead of match statements.
#[inline]
fn pack_bytes(bytes: &[u8]) -> u64 {
    bytes
        .iter()
        .fold(0u64, |acc, &b| (acc << 2) | PACK_TABLE[b as usize])
}

/// Unpacks packed bits to `Bytes` using a raw `usize` k.
///
/// Internal helper for use within the module where k is already validated
/// (e.g., from `bytes.len()` of a validated `Kmer`).
#[inline]
fn unpack_to_bytes_from_raw(packed_bits: u64, k: usize) -> Bytes {
    (0..k)
        .map(|i| {
            let shift = (k - 1 - i) * 2;
            let bits = ((packed_bits >> shift) & 0b11) as usize;
            UNPACK_TABLE[bits]
        })
        .collect()
}

// ============================================================================
// Zero-Copy Fast Path Functions
// ============================================================================

/// Computes the reverse complement of packed k-mer bits using pure arithmetic.
///
/// Each 2-bit pair is complemented (A<->T = 0<->3, C<->G = 1<->2, i.e. XOR with 3),
/// then the order of pairs is reversed.
///
/// # Arguments
///
/// * `bits` - Packed 2-bit representation of a k-mer
/// * `k` - K-mer length (number of bases, must be 1..=32)
///
/// # Examples
///
/// ```
/// use kmerust::kmer::reverse_complement_bits;
///
/// // ACGT (0b00_01_10_11) -> reverse complement ACGT (palindrome)
/// assert_eq!(reverse_complement_bits(0b00_01_10_11, 4), 0b00_01_10_11);
///
/// // AAA (0b00_00_00) -> TTT (0b11_11_11)
/// assert_eq!(reverse_complement_bits(0b00_00_00, 3), 0b11_11_11);
///
/// // TTT (0b11_11_11) -> AAA (0b00_00_00)
/// assert_eq!(reverse_complement_bits(0b11_11_11, 3), 0b00_00_00);
/// ```
///
/// # Panics
///
/// Panics in debug builds if `k` is 0 or greater than 32.
#[inline]
pub const fn reverse_complement_bits(bits: u64, k: usize) -> u64 {
    debug_assert!(k >= 1 && k <= 32, "k must be 1..=32");

    // Step 1: Complement all 2-bit pairs (A<->T, C<->G is XOR with 0b11 per pair)
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut rc = (!bits) & mask;

    // Step 2: Reverse the order of 2-bit pairs using parallel swap
    // Swap adjacent 2-bit groups
    rc = ((rc >> 2) & 0x3333_3333_3333_3333) | ((rc & 0x3333_3333_3333_3333) << 2);
    // Swap adjacent 4-bit groups
    rc = ((rc >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((rc & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    // Swap adjacent bytes
    rc = ((rc >> 8) & 0x00FF_00FF_00FF_00FF) | ((rc & 0x00FF_00FF_00FF_00FF) << 8);
    // Swap adjacent 2-byte groups
    rc = ((rc >> 16) & 0x0000_FFFF_0000_FFFF) | ((rc & 0x0000_FFFF_0000_FFFF) << 16);
    // Swap 4-byte halves
    rc = rc.rotate_left(32);

    // Step 3: Shift right to align k bases (they're now at the top of the u64)
    rc >> (64 - 2 * k)
}

/// Returns the canonical packed bits for a k-mer (min of forward and reverse complement).
///
/// This is a pure arithmetic operation with zero allocation.
#[inline]
pub const fn canonical_bits(packed: u64, k: usize) -> u64 {
    let rc = reverse_complement_bits(packed, k);
    if packed < rc {
        packed
    } else {
        rc
    }
}

/// Validates a DNA byte slice and packs it into a 64-bit integer in a single pass.
///
/// Returns the packed bits directly without any heap allocation.
/// Both uppercase and lowercase bases are accepted.
///
/// # Errors
///
/// Returns [`InvalidBaseError`] with the position and value of the first invalid byte.
///
/// # Examples
///
/// ```
/// use kmerust::kmer::validate_and_pack;
///
/// let bits = validate_and_pack(b"ACGT").unwrap();
/// assert_eq!(bits, 0b00_01_10_11);
///
/// // Lowercase works too
/// let bits = validate_and_pack(b"acgt").unwrap();
/// assert_eq!(bits, 0b00_01_10_11);
///
/// // Invalid bases are rejected
/// assert!(validate_and_pack(b"ACNT").is_err());
/// ```
#[inline]
pub fn validate_and_pack(seq: &[u8]) -> Result<u64, InvalidBaseError> {
    let mut packed: u64 = 0;
    for (i, &b) in seq.iter().enumerate() {
        match b {
            b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => {
                packed = (packed << 2) | PACK_TABLE[b as usize];
            }
            _ => {
                return Err(InvalidBaseError {
                    base: b,
                    position: i,
                });
            }
        }
    }
    Ok(packed)
}

/// Validates, packs, and canonicalizes a DNA byte slice in one fused operation.
///
/// This is the zero-allocation fast path for the k-mer counting hot loop.
/// Combines validation, 2-bit packing, and canonical selection without any
/// heap allocation.
///
/// # Errors
///
/// Returns [`InvalidBaseError`] with the position of the first invalid byte.
///
/// # Panics
///
/// Panics in debug builds if `seq` is empty.
///
/// # Examples
///
/// ```
/// use kmerust::kmer::pack_canonical;
///
/// // GATTACA and its RC TGTAATC -> GATTACA is smaller
/// let bits = pack_canonical(b"GATTACA").unwrap();
/// let rc_bits = pack_canonical(b"TGTAATC").unwrap();
/// assert_eq!(bits, rc_bits); // both map to same canonical form
/// ```
#[inline]
pub fn pack_canonical(seq: &[u8]) -> Result<u64, InvalidBaseError> {
    debug_assert!(!seq.is_empty(), "pack_canonical requires a non-empty slice");
    let packed = validate_and_pack(seq)?;
    Ok(canonical_bits(packed, seq.len()))
}

/// A single DNA base (nucleotide).
///
/// Used for converting between byte representations (ASCII) and
/// numeric representations (for bit packing).
pub enum KmerByte {
    /// Adenine
    A,
    /// Cytosine
    C,
    /// Guanine
    G,
    /// Thymine
    T,
}

impl From<&u8> for KmerByte {
    /// Converts a validated DNA base byte to `KmerByte`.
    ///
    /// # Safety Invariant
    ///
    /// This impl assumes the input byte has already been validated by
    /// `Kmer::from_sub`. Direct use with unvalidated input will panic in debug
    /// builds and produce undefined results in release builds.
    fn from(val: &u8) -> Self {
        debug_assert!(
            matches!(val, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'),
            "KmerByte::from called with invalid base: {val:#x}"
        );
        match val {
            b'A' | b'a' => Self::A,
            b'C' | b'c' => Self::C,
            b'G' | b'g' => Self::G,
            b'T' | b't' => Self::T,
            // SAFETY: Caller must ensure input is a valid DNA base.
            // This is guaranteed when called from Kmer methods after validation.
            _ => unreachable!("invalid base passed to KmerByte::from"),
        }
    }
}

impl From<KmerByte> for u8 {
    fn from(val: KmerByte) -> Self {
        match val {
            KmerByte::A => b'A',
            KmerByte::C => b'C',
            KmerByte::G => b'G',
            KmerByte::T => b'T',
        }
    }
}

impl From<u64> for KmerByte {
    /// Converts a 2-bit encoded value (0-3) to `KmerByte`.
    ///
    /// # Safety Invariant
    ///
    /// This impl assumes the input value is in range 0..=3, as produced by
    /// masking packed bits with `& 0b11`. Values outside this range will panic
    /// in debug builds.
    fn from(val: u64) -> Self {
        debug_assert!(
            val <= 3,
            "KmerByte::from called with invalid 2-bit value: {val}"
        );
        match val {
            0 => Self::A,
            1 => Self::C,
            2 => Self::G,
            3 => Self::T,
            // SAFETY: Caller must ensure input is a valid 2-bit value (0-3).
            // This is guaranteed when called from unpack functions with masked bits.
            _ => unreachable!("invalid 2-bit value passed to KmerByte::from"),
        }
    }
}

impl From<KmerByte> for u64 {
    fn from(val: KmerByte) -> Self {
        match val {
            KmerByte::A => 0,
            KmerByte::C => 1,
            KmerByte::G => 2,
            KmerByte::T => 3,
        }
    }
}

impl KmerByte {
    /// Fallibly converts a DNA base byte to `KmerByte`.
    ///
    /// Returns `Ok` for valid bases (A, C, G, T, case-insensitive),
    /// or `Err` for invalid bytes.
    ///
    /// # Example
    ///
    /// ```
    /// use kmerust::kmer::KmerByte;
    ///
    /// assert!(KmerByte::try_from_byte(b'A').is_ok());
    /// assert!(KmerByte::try_from_byte(b'N').is_err());
    /// ```
    pub const fn try_from_byte(val: u8) -> Result<Self, InvalidBaseError> {
        match val {
            b'A' | b'a' => Ok(Self::A),
            b'C' | b'c' => Ok(Self::C),
            b'G' | b'g' => Ok(Self::G),
            b'T' | b't' => Ok(Self::T),
            _ => Err(InvalidBaseError {
                base: val,
                position: 0,
            }),
        }
    }

    /// Fallibly converts a 2-bit encoded value to `KmerByte`.
    ///
    /// Returns `Ok` for values 0-3, or `Err` for values outside that range.
    ///
    /// # Example
    ///
    /// ```
    /// use kmerust::kmer::KmerByte;
    ///
    /// assert!(KmerByte::try_from_bits(0).is_ok()); // A
    /// assert!(KmerByte::try_from_bits(3).is_ok()); // T
    /// assert!(KmerByte::try_from_bits(4).is_err());
    /// ```
    #[allow(clippy::cast_possible_truncation)]
    pub const fn try_from_bits(val: u64) -> Result<Self, InvalidBaseError> {
        match val {
            0 => Ok(Self::A),
            1 => Ok(Self::C),
            2 => Ok(Self::G),
            3 => Ok(Self::T),
            // SAFETY: truncation is fine here as we're only interested in low bits
            _ => Err(InvalidBaseError {
                base: val as u8,
                position: 0,
            }),
        }
    }

    #[must_use]
    pub const fn reverse_complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
pub mod test {
    use super::*;

    #[test]
    fn bytes_from_valid_substring() {
        let sub = b"GATTACA";
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", kmer.bytes()), @r#"b"GATTACA""#);
    }

    #[test]
    fn from_substring_returns_err_for_invalid_substring() {
        let sub = b"N";
        let result = Kmer::from_sub(Bytes::copy_from_slice(sub));
        assert!(result.is_err());
    }

    #[test]
    fn from_sub_finds_invalid_byte_position() {
        let test_cases = [
            ("NACNN", 0, b'N'),
            ("ANCNG", 1, b'N'),
            ("AANTG", 2, b'N'),
            ("CCCNG", 3, b'N'),
            ("AACTN", 4, b'N'),
        ];

        for (dna, expected_pos, expected_base) in test_cases {
            let res = Kmer::from_sub(Bytes::copy_from_slice(dna.as_bytes()));
            let err = res.unwrap_err();
            assert_eq!(err.position, expected_pos, "for sequence {dna}");
            assert_eq!(err.base, expected_base, "for sequence {dna}");
        }
    }

    #[test]
    fn pack_unpack_roundtrip() {
        let sequences = ["ACGT", "AAAA", "TTTT", "CCCC", "GGGG", "GATTACA"];
        for seq in sequences {
            let kmer = Kmer::from_sub(Bytes::copy_from_slice(seq.as_bytes())).unwrap();
            let packed = kmer.pack();
            let k = KmerLength::new(seq.len()).unwrap();
            let unpacked = unpack_to_bytes(packed.packed_bits(), k);
            assert_eq!(unpacked.as_ref(), seq.as_bytes());
        }
    }

    #[test]
    fn pack_unpack_roundtrip_various_lengths() {
        for k_val in 1..=32 {
            let seq = "A".repeat(k_val);
            let kmer = Kmer::from_sub(Bytes::copy_from_slice(seq.as_bytes())).unwrap();
            let packed = kmer.pack();
            let k = KmerLength::new(k_val).unwrap();
            let unpacked = unpack_to_bytes(packed.packed_bits(), k);
            assert_eq!(unpacked.as_ref(), seq.as_bytes());
        }
    }

    #[test]
    fn canonical_selects_lexicographically_smaller() {
        // ACGT and ACGT (reverse complement) - ACGT is palindromic
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(b"ACGT"))
            .unwrap()
            .pack()
            .canonical();
        assert_eq!(kmer.bytes().as_ref(), b"ACGT");
        assert!(!kmer.is_reverse_complement());

        // AAA -> TTT reverse complement, AAA is smaller
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(b"AAA"))
            .unwrap()
            .pack()
            .canonical();
        assert_eq!(kmer.bytes().as_ref(), b"AAA");
        assert!(!kmer.is_reverse_complement());

        // TTT -> AAA reverse complement, AAA is smaller
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(b"TTT"))
            .unwrap()
            .pack()
            .canonical();
        assert_eq!(kmer.bytes().as_ref(), b"AAA");
        assert!(kmer.is_reverse_complement());

        // GATTACA -> TGTAATC reverse complement, GATTACA is smaller
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(b"GATTACA"))
            .unwrap()
            .pack()
            .canonical();
        assert_eq!(kmer.bytes().as_ref(), b"GATTACA");
        assert!(!kmer.is_reverse_complement());

        // TGTAATC -> GATTACA reverse complement, GATTACA is smaller
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(b"TGTAATC"))
            .unwrap()
            .pack()
            .canonical();
        assert_eq!(kmer.bytes().as_ref(), b"GATTACA");
        assert!(kmer.is_reverse_complement());
    }

    #[test]
    fn kmer_byte_reverse_complement() {
        assert!(matches!(KmerByte::A.reverse_complement(), KmerByte::T));
        assert!(matches!(KmerByte::T.reverse_complement(), KmerByte::A));
        assert!(matches!(KmerByte::C.reverse_complement(), KmerByte::G));
        assert!(matches!(KmerByte::G.reverse_complement(), KmerByte::C));
    }

    #[test]
    fn kmer_byte_to_u64_roundtrip() {
        for (byte, expected) in [(b'A', 0u64), (b'C', 1u64), (b'G', 2u64), (b'T', 3u64)] {
            let kmer_byte: KmerByte = (&byte).into();
            let val: u64 = kmer_byte.into();
            assert_eq!(val, expected);
            let back: KmerByte = val.into();
            let back_byte: u8 = back.into();
            assert_eq!(back_byte, byte);
        }
    }

    #[test]
    fn empty_kmer() {
        let kmer = Kmer::from_sub(Bytes::new()).unwrap();
        assert!(kmer.bytes().is_empty());
    }

    #[test]
    fn single_base_kmers() {
        for base in [b'A', b'C', b'G', b'T'] {
            let kmer = Kmer::from_sub(Bytes::copy_from_slice(&[base])).unwrap();
            assert_eq!(kmer.bytes().as_ref(), &[base]);
        }
    }

    #[test]
    fn soft_masked_bases_converted_to_uppercase() {
        // Test that lowercase bases (soft-masked) are accepted and converted
        let sub = b"AAAa";
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", kmer.bytes()), @r#"b"AAAA""#);
    }

    #[test]
    fn soft_masked_all_lowercase() {
        // Test all lowercase sequence
        let sub = b"gattaca";
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", kmer.bytes()), @r#"b"GATTACA""#);
    }

    #[test]
    fn soft_masked_mixed_case() {
        // Test mixed case sequence
        let sub = b"AcGt";
        let kmer = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", kmer.bytes()), @r#"b"ACGT""#);
    }

    #[test]
    fn kmer_length_valid_range() {
        for k in 1..=32 {
            let result = KmerLength::new(k);
            assert!(result.is_ok(), "k={k} should be valid");
            assert_eq!(result.unwrap().get(), k);
        }
    }

    #[test]
    fn kmer_length_rejects_zero() {
        let result = KmerLength::new(0);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert_eq!(err.k, 0);
        assert_eq!(err.min, KmerLength::MIN);
        assert_eq!(err.max, KmerLength::MAX);
    }

    #[test]
    fn kmer_length_rejects_too_large() {
        for k in [33, 64, 100, 1000] {
            let result = KmerLength::new(k);
            assert!(result.is_err(), "k={k} should be invalid");
        }
    }

    #[test]
    fn kmer_length_try_from() {
        let k: Result<KmerLength, _> = 21usize.try_into();
        assert!(k.is_ok());
        assert_eq!(k.unwrap().get(), 21);
    }

    #[test]
    fn kmer_length_into_usize() {
        let k = KmerLength::new(21).unwrap();
        let n: usize = k.into();
        assert_eq!(n, 21);
    }

    #[test]
    fn kmer_length_as_u8() {
        let k = KmerLength::new(21).unwrap();
        assert_eq!(k.as_u8(), 21);
    }

    #[test]
    fn unpack_to_string_works() {
        let k = KmerLength::new(4).unwrap();
        // ACGT = 00 01 10 11 = 27
        let s = unpack_to_string(0b00_01_10_11, k);
        assert_eq!(s, "ACGT");
    }

    #[test]
    fn reverse_complement_bits_single_base() {
        // A(00) -> T(11), T(11) -> A(00), C(01) -> G(10), G(10) -> C(01)
        assert_eq!(reverse_complement_bits(0b00, 1), 0b11);
        assert_eq!(reverse_complement_bits(0b11, 1), 0b00);
        assert_eq!(reverse_complement_bits(0b01, 1), 0b10);
        assert_eq!(reverse_complement_bits(0b10, 1), 0b01);
    }

    #[test]
    fn reverse_complement_bits_max_k() {
        // k=32: all A's -> all T's
        let all_a = 0u64;
        let all_t = u64::MAX; // 32 T bases = all 1s
        assert_eq!(reverse_complement_bits(all_a, 32), all_t);
        assert_eq!(reverse_complement_bits(all_t, 32), all_a);
    }

    #[test]
    fn reverse_complement_bits_palindrome() {
        // ACGT is a palindrome (RC = ACGT)
        assert_eq!(reverse_complement_bits(0b00_01_10_11, 4), 0b00_01_10_11);

        // AT is a palindrome (RC = AT)
        assert_eq!(reverse_complement_bits(0b00_11, 2), 0b00_11);
    }

    #[test]
    fn canonical_bits_returns_smaller() {
        // AAA(0b000000) vs TTT(0b111111) -> AAA wins
        assert_eq!(canonical_bits(0b00_00_00, 3), 0b00_00_00);
        assert_eq!(canonical_bits(0b11_11_11, 3), 0b00_00_00);
    }

    #[test]
    fn canonical_bits_palindrome_returns_forward() {
        // ACGT is palindromic, forward == RC, should return forward value
        let acgt = 0b00_01_10_11u64;
        assert_eq!(canonical_bits(acgt, 4), acgt);
    }

    #[test]
    #[should_panic(expected = "pack_canonical requires a non-empty slice")]
    fn pack_canonical_panics_on_empty_slice() {
        let _ = pack_canonical(b"");
    }

    #[test]
    fn validate_and_pack_rejects_invalid() {
        let err = validate_and_pack(b"ACN").unwrap_err();
        assert_eq!(err.position, 2);
        assert_eq!(err.base, b'N');
    }

    #[test]
    fn pack_canonical_matches_type_state_pipeline() {
        let sequences = [
            "GATTACA",
            "TGTAATC",
            "AAA",
            "TTT",
            "ACGT",
            "A",
            "T",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // k=32
        ];
        for seq in sequences {
            let via_pipeline = Kmer::from_sub(Bytes::copy_from_slice(seq.as_bytes()))
                .unwrap()
                .pack()
                .canonical()
                .packed_bits();
            let via_fast_path = pack_canonical(seq.as_bytes()).unwrap();
            assert_eq!(via_pipeline, via_fast_path, "mismatch for {seq}");
        }
    }

    #[test]
    fn type_state_flow() {
        // Demonstrate the type-state flow
        let unpacked: Kmer<Unpacked> = Kmer::from_sub(Bytes::from_static(b"GATTACA")).unwrap();
        let packed: Kmer<Packed> = unpacked.pack();
        let canonical: Kmer<Canonical> = packed.canonical();

        // Can access packed_bits on both Packed and Canonical
        assert!(canonical.packed_bits() > 0);
        assert!(!canonical.is_reverse_complement());
    }
}
