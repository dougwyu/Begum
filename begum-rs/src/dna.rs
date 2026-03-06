/// Returns true if `base` is compatible with the IUPAC ambiguity code `iupac`.
/// Both inputs are treated case-insensitively.
pub fn iupac_matches(iupac: u8, base: u8) -> bool {
    let base = base.to_ascii_uppercase();
    match iupac.to_ascii_uppercase() {
        b'A' => base == b'A',
        b'C' => base == b'C',
        b'G' => base == b'G',
        b'T' | b'U' => base == b'T' || base == b'U',
        b'R' => base == b'A' || base == b'G',
        b'Y' => base == b'C' || base == b'T',
        b'S' => base == b'G' || base == b'C',
        b'W' => base == b'A' || base == b'T',
        b'K' => base == b'G' || base == b'T',
        b'M' => base == b'A' || base == b'C',
        b'B' => base == b'C' || base == b'G' || base == b'T',
        b'D' => base == b'A' || base == b'G' || base == b'T',
        b'H' => base == b'A' || base == b'C' || base == b'T',
        b'V' => base == b'A' || base == b'C' || base == b'G',
        b'N' | b'I' => true,
        _ => false,
    }
}

/// Count positions where primer[i] does not match seq_slice[i] under IUPAC rules.
/// primer and seq_slice must have the same length.
pub fn count_mismatches(primer: &[u8], seq_slice: &[u8]) -> usize {
    primer
        .iter()
        .zip(seq_slice.iter())
        .filter(|(&p, &s)| !iupac_matches(p, s))
        .count()
}

/// Return all (start, end) positions where primer occurs in seq
/// with at most max_mismatches substitutions (IUPAC-aware).
pub fn find_primer_occurrences(
    primer: &[u8],
    seq: &[u8],
    max_mismatches: usize,
) -> Vec<(usize, usize)> {
    let plen = primer.len();
    let slen = seq.len();
    if plen > slen {
        return vec![];
    }
    (0..=(slen - plen))
        .filter(|&i| count_mismatches(primer, &seq[i..i + plen]) <= max_mismatches)
        .map(|i| (i, i + plen))
        .collect()
}

/// Return (start, end, total_count) of the FIRST occurrence of primer in seq.
/// Returns (-1, -1, 0) if no match. total_count is the number of all matches found.
/// Mirrors Python dna_utility.find_first_match.
pub fn find_primer_first(
    primer: &[u8],
    seq: &[u8],
    max_mismatches: usize,
) -> (i64, i64, usize) {
    let hits = find_primer_occurrences(primer, seq, max_mismatches);
    match hits.first() {
        Some(&(s, e)) => (s as i64, e as i64, hits.len()),
        None => (-1, -1, 0),
    }
}

/// Return (start, end, total_count) of the LAST occurrence of primer in seq.
/// Returns (-1, -1, 0) if no match. total_count is the number of all matches found.
/// Mirrors Python dna_utility.find_last_match.
pub fn find_primer_last(
    primer: &[u8],
    seq: &[u8],
    max_mismatches: usize,
) -> (i64, i64, usize) {
    let hits = find_primer_occurrences(primer, seq, max_mismatches);
    match hits.last() {
        Some(&(s, e)) => (s as i64, e as i64, hits.len()),
        None => (-1, -1, 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iupac_exact() {
        assert!(iupac_matches(b'A', b'A'));
        assert!(!iupac_matches(b'A', b'C'));
        assert!(iupac_matches(b'G', b'G'));
        assert!(!iupac_matches(b'T', b'G'));
    }

    #[test]
    fn test_iupac_ambiguous() {
        assert!(iupac_matches(b'R', b'A')); // R = A or G
        assert!(iupac_matches(b'R', b'G'));
        assert!(!iupac_matches(b'R', b'C'));
        assert!(iupac_matches(b'Y', b'C')); // Y = C or T
        assert!(iupac_matches(b'Y', b'T'));
        assert!(!iupac_matches(b'Y', b'A'));
        assert!(iupac_matches(b'N', b'T')); // N = any
        assert!(iupac_matches(b'N', b'A'));
        assert!(iupac_matches(b'I', b'G')); // I = any (like N)
    }

    #[test]
    fn test_iupac_case_insensitive() {
        assert!(iupac_matches(b'r', b'a')); // lowercase works too
        assert!(iupac_matches(b'R', b'a'));
        assert!(iupac_matches(b'r', b'A'));
    }

    #[test]
    fn test_count_mismatches_exact() {
        assert_eq!(count_mismatches(b"GCATGC", b"GCATGC"), 0);
        assert_eq!(count_mismatches(b"GCATGC", b"GCATGX"), 1);
        assert_eq!(count_mismatches(b"GCATGC", b"XXXXXX"), 6);
    }

    #[test]
    fn test_count_mismatches_iupac() {
        // R matches A or G — so RCAT vs ACAT = 0 mismatches
        assert_eq!(count_mismatches(b"RCATGC", b"ACATGC"), 0);
        assert_eq!(count_mismatches(b"RCATGC", b"GCATGC"), 0);
        assert_eq!(count_mismatches(b"RCATGC", b"CCATGC"), 1); // C not in R
    }

    #[test]
    fn test_find_primer_first_exact() {
        // primer at position 7 in: AACCGGT[GCATGC]AAATTT
        let seq = b"AACCGGTGCATGCAAATTT";
        let primer = b"GCATGC";
        let (start, end, count) = find_primer_first(primer, seq, 0);
        assert_eq!(start, 7);
        assert_eq!(end, 13);
        assert_eq!(count, 1);
    }

    #[test]
    fn test_find_primer_first_no_match() {
        let (start, _end, count) = find_primer_first(b"TTTTTT", b"AACCGGT", 0);
        assert_eq!(start, -1);
        assert_eq!(count, 0);
    }

    #[test]
    fn test_find_primer_first_with_mismatch() {
        // 1 mismatch allowed — GCATGX matches GCATGC with 1 sub
        let seq = b"AAAGCATGXBBB";
        let primer = b"GCATGC";
        let (start, end, count) = find_primer_first(primer, seq, 1);
        assert_eq!(start, 3);
        assert_eq!(end, 9);
        assert_eq!(count, 1);
    }

    #[test]
    fn test_find_primer_last() {
        // primer appears twice — last match should be the second one
        let seq = b"GCATGCAAAGCATGCBBB";
        let primer = b"GCATGC";
        let (start, end, count) = find_primer_last(primer, seq, 0);
        assert_eq!(start, 9);
        assert_eq!(end, 15);
        assert_eq!(count, 2);
    }

    #[test]
    fn test_find_primer_last_no_match() {
        let (start, _end, count) = find_primer_last(b"TTTTTT", b"AACCGGT", 0);
        assert_eq!(start, -1);
        assert_eq!(count, 0);
    }

    #[test]
    fn test_primer_longer_than_seq() {
        let (start, _end, count) = find_primer_first(b"GCATGCGCATGC", b"GCATGC", 0);
        assert_eq!(start, -1);
        assert_eq!(count, 0);
    }
}
