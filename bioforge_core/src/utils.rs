//! Utility functions for BioForge Core.

use std::path::Path;

/// Check if file is gzipped.
pub fn is_gzipped<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref()
        .extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false)
}

/// Get file extension without .gz.
pub fn get_base_extension<P: AsRef<Path>>(path: P) -> Option<String> {
    let path = path.as_ref();
    let ext = path.extension()?.to_str()?;
    
    if ext == "gz" {
        // Get the extension before .gz
        let stem = path.file_stem()?;
        Path::new(stem)
            .extension()
            .and_then(|e| e.to_str())
            .map(|s| s.to_string())
    } else {
        Some(ext.to_string())
    }
}

/// Reverse complement a DNA sequence.
pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            'N' | 'n' => 'N',
            _ => 'N',
        })
        .collect()
}

/// Calculate GC content of a sequence.
pub fn gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq.iter().filter(|b| **b == b'G' || **b == b'C').count();
    gc as f64 / seq.len() as f64 * 100.0
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement(""), "");
    }
    
    #[test]
    fn test_gc_content() {
        assert_eq!(gc_content(b"GCGC"), 100.0);
        assert_eq!(gc_content(b"ATAT"), 0.0);
        assert_eq!(gc_content(b"ATGC"), 50.0);
    }
    
    #[test]
    fn test_is_gzipped() {
        assert!(is_gzipped("file.fastq.gz"));
        assert!(!is_gzipped("file.fastq"));
    }
}
