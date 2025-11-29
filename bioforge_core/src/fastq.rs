//! FASTQ parsing and QC module.
//!
//! High-performance FASTQ parsing using needletail.

use pyo3::prelude::*;
use needletail::parse_fastx_file;

/// FASTQ record representation.
#[pyclass]
#[derive(Clone)]
pub struct FastqRecord {
    #[pyo3(get)]
    pub id: String,
    #[pyo3(get)]
    pub sequence: String,
    #[pyo3(get)]
    pub quality: String,
}

#[pymethods]
impl FastqRecord {
    #[new]
    fn new(id: String, sequence: String, quality: String) -> Self {
        FastqRecord { id, sequence, quality }
    }
    
    /// Get sequence length.
    fn len(&self) -> usize {
        self.sequence.len()
    }
    
    /// Calculate mean quality score.
    fn mean_quality(&self) -> f64 {
        if self.quality.is_empty() {
            return 0.0;
        }
        let sum: u64 = self.quality.bytes().map(|b| (b - 33) as u64).sum();
        sum as f64 / self.quality.len() as f64
    }
    
    /// Calculate GC content.
    fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc_count = self.sequence.chars().filter(|c| *c == 'G' || *c == 'C').count();
        gc_count as f64 / self.sequence.len() as f64 * 100.0
    }
}

/// FASTQ statistics.
#[pyclass]
#[derive(Clone, Default)]
pub struct FastqStats {
    #[pyo3(get)]
    pub total_reads: u64,
    #[pyo3(get)]
    pub total_bases: u64,
    #[pyo3(get)]
    pub mean_length: f64,
    #[pyo3(get)]
    pub min_length: usize,
    #[pyo3(get)]
    pub max_length: usize,
    #[pyo3(get)]
    pub mean_quality: f64,
    #[pyo3(get)]
    pub gc_content: f64,
    #[pyo3(get)]
    pub q20_bases: u64,
    #[pyo3(get)]
    pub q30_bases: u64,
}

#[pymethods]
impl FastqStats {
    fn __repr__(&self) -> String {
        format!(
            "FastqStats(reads={}, bases={}, mean_len={:.1}, gc={:.1}%, q30={:.1}%)",
            self.total_reads,
            self.total_bases,
            self.mean_length,
            self.gc_content,
            if self.total_bases > 0 { self.q30_bases as f64 / self.total_bases as f64 * 100.0 } else { 0.0 }
        )
    }
}

/// FASTQ file reader.
#[pyclass]
pub struct FastqReader {
    path: String,
}

#[pymethods]
impl FastqReader {
    #[new]
    fn new(path: String) -> PyResult<Self> {
        Ok(FastqReader { path })
    }
    
    /// Get statistics for the file.
    fn stats(&self) -> PyResult<FastqStats> {
        calculate_stats(&self.path)
    }
    
    /// Read all records into a list.
    fn read_all(&self) -> PyResult<Vec<FastqRecord>> {
        let mut records = Vec::new();
        let mut reader = parse_fastx_file(&self.path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        
        while let Some(Ok(record)) = reader.next() {
            let id = String::from_utf8_lossy(record.id()).to_string();
            let seq = record.seq();
            let sequence = String::from_utf8_lossy(&seq).to_string();
            let quality = record.qual()
                .map(|q| String::from_utf8_lossy(q).to_string())
                .unwrap_or_default();
            records.push(FastqRecord { id, sequence, quality });
        }
        
        Ok(records)
    }
}

/// Count reads in a FASTQ file.
#[pyfunction]
pub fn count_reads(path: &str) -> PyResult<u64> {
    let mut reader = parse_fastx_file(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    
    let mut count = 0u64;
    while let Some(Ok(_)) = reader.next() {
        count += 1;
    }
    Ok(count)
}

/// Calculate statistics for a FASTQ file.
#[pyfunction]
pub fn calculate_stats(path: &str) -> PyResult<FastqStats> {
    let mut reader = parse_fastx_file(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    
    let mut stats = FastqStats::default();
    stats.min_length = usize::MAX;
    
    let mut total_quality_sum: u64 = 0;
    let mut gc_count: u64 = 0;
    
    while let Some(Ok(record)) = reader.next() {
        let seq = record.seq();
        let len = seq.len();
        
        stats.total_reads += 1;
        stats.total_bases += len as u64;
        
        if len < stats.min_length {
            stats.min_length = len;
        }
        if len > stats.max_length {
            stats.max_length = len;
        }
        
        // GC content
        gc_count += seq.iter().filter(|b| **b == b'G' || **b == b'C').count() as u64;
        
        // Quality scores
        if let Some(qual) = record.qual() {
            for &q in qual {
                let score = q - 33;
                total_quality_sum += score as u64;
                if score >= 20 {
                    stats.q20_bases += 1;
                }
                if score >= 30 {
                    stats.q30_bases += 1;
                }
            }
        }
    }
    
    if stats.total_reads > 0 {
        stats.mean_length = stats.total_bases as f64 / stats.total_reads as f64;
        stats.gc_content = gc_count as f64 / stats.total_bases as f64 * 100.0;
        stats.mean_quality = total_quality_sum as f64 / stats.total_bases as f64;
    }
    
    if stats.min_length == usize::MAX {
        stats.min_length = 0;
    }
    
    Ok(stats)
}
