//! BAM/SAM alignment operations module.
//!
//! High-performance BAM reading and statistics using noodles.

use pyo3::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::bam;

/// Alignment statistics.
#[pyclass]
#[derive(Clone, Default, Debug)]
pub struct AlignmentStats {
    #[pyo3(get)]
    pub total_reads: u64,
    #[pyo3(get)]
    pub mapped_reads: u64,
    #[pyo3(get)]
    pub unmapped_reads: u64,
    #[pyo3(get)]
    pub properly_paired: u64,
    #[pyo3(get)]
    pub duplicates: u64,
    #[pyo3(get)]
    pub secondary: u64,
    #[pyo3(get)]
    pub supplementary: u64,
    #[pyo3(get)]
    pub mapping_rate: f64,
    #[pyo3(get)]
    pub mean_mapq: f64,
    #[pyo3(get)]
    pub mean_insert_size: f64,
    #[pyo3(get)]
    pub reads_per_chromosome: HashMap<String, u64>,
}

#[pymethods]
impl AlignmentStats {
    #[new]
    fn new() -> Self {
        AlignmentStats::default()
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentStats(total={}, mapped={}, rate={:.2}%, mapq={:.1})",
            self.total_reads, self.mapped_reads, self.mapping_rate, self.mean_mapq
        )
    }

    /// Get summary as dict
    fn to_dict(&self) -> HashMap<String, f64> {
        let mut d = HashMap::new();
        d.insert("total_reads".into(), self.total_reads as f64);
        d.insert("mapped_reads".into(), self.mapped_reads as f64);
        d.insert("unmapped_reads".into(), self.unmapped_reads as f64);
        d.insert("properly_paired".into(), self.properly_paired as f64);
        d.insert("duplicates".into(), self.duplicates as f64);
        d.insert("mapping_rate".into(), self.mapping_rate);
        d.insert("mean_mapq".into(), self.mean_mapq);
        d
    }
}

/// BAM file reader with statistics and region queries.
#[pyclass]
pub struct BamReader {
    path: String,
}

#[pymethods]
impl BamReader {
    #[new]
    fn new(path: String) -> PyResult<Self> {
        // Validate file exists
        if !Path::new(&path).exists() {
            return Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
                format!("BAM file not found: {}", path),
            ));
        }

        Ok(BamReader { path })
    }

    /// Get alignment statistics for the BAM file.
    fn stats(&self) -> PyResult<AlignmentStats> {
        calculate_bam_stats(&self.path)
    }

    /// Get coverage for a genomic region.
    fn coverage(&self, chrom: &str, start: u64, end: u64) -> PyResult<Vec<u32>> {
        get_coverage(&self.path, chrom, start, end)
    }

    /// Count reads per chromosome.
    fn reads_per_chromosome(&self) -> PyResult<HashMap<String, u64>> {
        let stats = calculate_bam_stats(&self.path)?;
        Ok(stats.reads_per_chromosome)
    }

    /// Get reference sequences from header.
    fn reference_sequences(&self) -> PyResult<Vec<(String, u64)>> {
        let file = File::open(&self.path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let mut reader = bam::io::Reader::new(BufReader::new(file));

        let header = reader
            .read_header()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        let refs: Vec<(String, u64)> = header
            .reference_sequences()
            .iter()
            .map(|(name, seq)| {
                let len = seq.length().get() as u64;
                (name.to_string(), len)
            })
            .collect();

        Ok(refs)
    }

    /// Get the file path.
    fn get_path(&self) -> &str {
        &self.path
    }

    /// Check if BAM has index.
    fn has_index(&self) -> bool {
        let bai_path = format!("{}.bai", self.path);
        Path::new(&bai_path).exists()
    }
}

/// Calculate comprehensive BAM statistics.
fn calculate_bam_stats(bam_path: &str) -> PyResult<AlignmentStats> {
    let file = File::open(bam_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    let mut reader = bam::io::Reader::new(BufReader::new(file));

    let header = reader
        .read_header()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut stats = AlignmentStats::default();
    let mut total_mapq: u64 = 0;
    let mut total_insert_size: i64 = 0;
    let mut insert_count: u64 = 0;
    let mut reads_per_chr: HashMap<String, u64> = HashMap::new();

    // Build reference name lookup
    let ref_names: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| name.to_string())
        .collect();

    for result in reader.records() {
        let record = result
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        stats.total_reads += 1;

        // Get flags
        let flags = record.flags();

        // Check mapping status
        if !flags.is_unmapped() {
            stats.mapped_reads += 1;

            // Get MAPQ
            if let Some(mapq) = record.mapping_quality() {
                total_mapq += mapq.get() as u64;
            }

            // Count reads per chromosome
            if let Some(ref_id_result) = record.reference_sequence_id() {
                if let Ok(ref_id) = ref_id_result {
                    if let Some(chr_name) = ref_names.get(ref_id) {
                        *reads_per_chr.entry(chr_name.clone()).or_insert(0) += 1;
                    }
                }
            }
        } else {
            stats.unmapped_reads += 1;
        }

        // Check pair status
        if flags.is_properly_segmented() {
            stats.properly_paired += 1;
        }

        // Check duplicate
        if flags.is_duplicate() {
            stats.duplicates += 1;
        }

        // Check secondary/supplementary
        if flags.is_secondary() {
            stats.secondary += 1;
        }
        if flags.is_supplementary() {
            stats.supplementary += 1;
        }

        // Insert size (for properly paired reads)
        if flags.is_properly_segmented() && !flags.is_unmapped() {
            let tlen = record.template_length();
            if tlen > 0 {
                total_insert_size += tlen as i64;
                insert_count += 1;
            }
        }
    }

    // Calculate rates
    if stats.total_reads > 0 {
        stats.mapping_rate = (stats.mapped_reads as f64 / stats.total_reads as f64) * 100.0;
    }

    if stats.mapped_reads > 0 {
        stats.mean_mapq = total_mapq as f64 / stats.mapped_reads as f64;
    }

    if insert_count > 0 {
        stats.mean_insert_size = total_insert_size as f64 / insert_count as f64;
    }

    stats.reads_per_chromosome = reads_per_chr;

    Ok(stats)
}

/// Get coverage for a genomic region.
#[pyfunction]
pub fn get_coverage(bam_path: &str, chrom: &str, start: u64, end: u64) -> PyResult<Vec<u32>> {
    let file = File::open(bam_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    let mut reader = bam::io::Reader::new(BufReader::new(file));

    let header = reader
        .read_header()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    // Find reference sequence ID
    let ref_id = header
        .reference_sequences()
        .iter()
        .position(|(name, _)| name.to_string() == chrom);

    let ref_id = match ref_id {
        Some(id) => id,
        None => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Chromosome not found: {}",
                chrom
            )))
        }
    };

    let region_len = (end - start) as usize;
    let mut coverage = vec![0u32; region_len];

    for result in reader.records() {
        let record = result
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        // Skip unmapped reads
        if record.flags().is_unmapped() {
            continue;
        }

        // Check if on target chromosome
        if let Some(rec_ref_id_result) = record.reference_sequence_id() {
            if let Ok(rec_ref_id) = rec_ref_id_result {
                if rec_ref_id != ref_id {
                    continue;
                }
            } else {
                continue;
            }
        } else {
            continue;
        }

        // Get alignment position
        let aln_start_result = record.alignment_start();
        let aln_start = match aln_start_result {
            Some(Ok(pos)) => pos.get() as u64,
            _ => continue,
        };

        // Get CIGAR and calculate covered positions
        let cigar = record.cigar();
        let mut ref_pos = aln_start;

        for op_result in cigar.iter() {
            let op = match op_result {
                Ok(o) => o,
                Err(_) => continue,
            };

            use noodles::sam::alignment::record::cigar::op::Kind;
            let len = op.len();

            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for i in 0..len {
                        let pos = ref_pos + i as u64;
                        if pos >= start && pos < end {
                            let idx = (pos - start) as usize;
                            coverage[idx] += 1;
                        }
                    }
                    ref_pos += len as u64;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len as u64;
                }
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {
                    // These don't consume reference
                }
            }
        }
    }

    Ok(coverage)
}

/// Calculate mean coverage for a BAM file.
#[pyfunction]
pub fn mean_coverage(bam_path: &str) -> PyResult<f64> {
    let stats = calculate_bam_stats(bam_path)?;

    // Simple estimation: mapped_reads * avg_read_length / genome_size
    // This is a rough estimate without full depth calculation
    Ok(stats.mapped_reads as f64 * 150.0 / 3_000_000_000.0)
}

/// Count mapped reads in a BAM file.
#[pyfunction]
pub fn count_mapped_reads(bam_path: &str) -> PyResult<u64> {
    let stats = calculate_bam_stats(bam_path)?;
    Ok(stats.mapped_reads)
}
