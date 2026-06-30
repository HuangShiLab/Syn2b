//! Utility functions for 2bSyn

/// Reverse complement of a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => b'N',
        })
        .collect()
}

/// Check if a byte sequence contains only valid DNA bases (ACGTN)
pub fn is_valid_dna(seq: &[u8]) -> bool {
    seq.iter().all(|&b| matches!(b, b'A'|b'C'|b'G'|b'T'|b'a'|b'c'|b'g'|b't'|b'N'|b'n'))
}

/// GC content of a sequence (0.0 to 1.0)
pub fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() { return 0.0; }
    let gc = seq.iter().filter(|&&b| matches!(b, b'G'|b'g'|b'C'|b'c')).count();
    gc as f64 / seq.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATCG"), vec![b'C', b'G', b'A', b'T']);
        assert_eq!(reverse_complement(b"AATT"), vec![b'A', b'A', b'T', b'T']);
    }

    #[test]
    fn test_is_valid_dna() {
        assert!(is_valid_dna(b"ATCG"));
        assert!(!is_valid_dna(b"XYZ"));
    }

    #[test]
    fn test_gc_content() {
        assert!((gc_content(b"GCGC") - 1.0).abs() < 1e-9);
        assert!((gc_content(b"ATAT") - 0.0).abs() < 1e-9);
        assert!((gc_content(b"ATCG") - 0.5).abs() < 1e-9);
    }
}
