//! VCF parsing and variant operations module.
//!
//! High-performance VCF operations using noodles.

use pyo3::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

/// Variant representation.
#[pyclass]
#[derive(Clone, Debug)]
pub struct Variant {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub pos: u64,
    #[pyo3(get)]
    pub id: Option<String>,
    #[pyo3(get)]
    pub ref_allele: String,
    #[pyo3(get)]
    pub alt_allele: String,
    #[pyo3(get)]
    pub qual: Option<f64>,
    #[pyo3(get)]
    pub filter: String,
    #[pyo3(get)]
    pub info: HashMap<String, String>,
    #[pyo3(get)]
    pub genotype: Option<String>,
}

#[pymethods]
impl Variant {
    #[new]
    #[pyo3(signature = (chrom, pos, ref_allele, alt_allele, id=None, qual=None, filter=None))]
    fn new(
        chrom: String,
        pos: u64,
        ref_allele: String,
        alt_allele: String,
        id: Option<String>,
        qual: Option<f64>,
        filter: Option<String>,
    ) -> Self {
        Variant {
            chrom,
            pos,
            id,
            ref_allele,
            alt_allele,
            qual,
            filter: filter.unwrap_or_else(|| "PASS".to_string()),
            info: HashMap::new(),
            genotype: None,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Variant({}:{} {}>{}, filter={}, qual={:?})",
            self.chrom, self.pos, self.ref_allele, self.alt_allele, self.filter, self.qual
        )
    }

    /// Check if variant is a SNP (single nucleotide polymorphism).
    fn is_snp(&self) -> bool {
        self.ref_allele.len() == 1 && self.alt_allele.len() == 1
    }

    /// Check if variant is an indel (insertion or deletion).
    fn is_indel(&self) -> bool {
        self.ref_allele.len() != self.alt_allele.len()
    }

    /// Check if variant is an insertion.
    fn is_insertion(&self) -> bool {
        self.ref_allele.len() < self.alt_allele.len()
    }

    /// Check if variant is a deletion.
    fn is_deletion(&self) -> bool {
        self.ref_allele.len() > self.alt_allele.len()
    }

    /// Get variant type string.
    fn variant_type(&self) -> String {
        if self.is_snp() {
            "SNV".to_string()
        } else if self.is_insertion() {
            "INS".to_string()
        } else if self.is_deletion() {
            "DEL".to_string()
        } else {
            "MNV".to_string()
        }
    }

    /// Get indel length (positive for insertions, negative for deletions).
    fn indel_length(&self) -> i32 {
        self.alt_allele.len() as i32 - self.ref_allele.len() as i32
    }

    /// Check if variant passes filters.
    fn is_pass(&self) -> bool {
        self.filter == "PASS" || self.filter == "."
    }

    /// Check if genotype is homozygous.
    fn is_homozygous(&self) -> bool {
        if let Some(ref gt) = self.genotype {
            let parts: Vec<&str> = gt.split(['/', '|'].as_ref()).collect();
            parts.len() >= 2 && parts[0] == parts[1] && parts[0] != "0"
        } else {
            false
        }
    }

    /// Check if genotype is heterozygous.
    fn is_heterozygous(&self) -> bool {
        if let Some(ref gt) = self.genotype {
            let parts: Vec<&str> = gt.split(['/', '|'].as_ref()).collect();
            parts.len() >= 2 && parts[0] != parts[1]
        } else {
            false
        }
    }

    /// Get info field value.
    fn get_info(&self, key: &str) -> Option<String> {
        self.info.get(key).cloned()
    }
}

/// VCF statistics.
#[pyclass]
#[derive(Clone, Default, Debug)]
pub struct VcfStats {
    #[pyo3(get)]
    pub total_variants: u64,
    #[pyo3(get)]
    pub snps: u64,
    #[pyo3(get)]
    pub indels: u64,
    #[pyo3(get)]
    pub insertions: u64,
    #[pyo3(get)]
    pub deletions: u64,
    #[pyo3(get)]
    pub multi_allelic: u64,
    #[pyo3(get)]
    pub pass_filter: u64,
    #[pyo3(get)]
    pub filtered: u64,
    #[pyo3(get)]
    pub transitions: u64,
    #[pyo3(get)]
    pub transversions: u64,
    #[pyo3(get)]
    pub ti_tv_ratio: f64,
    #[pyo3(get)]
    pub het_hom_ratio: f64,
    #[pyo3(get)]
    pub mean_quality: f64,
}

#[pymethods]
impl VcfStats {
    fn __repr__(&self) -> String {
        format!(
            "VcfStats(total={}, snps={}, indels={}, ti_tv={:.2})",
            self.total_variants, self.snps, self.indels, self.ti_tv_ratio
        )
    }

    fn to_dict(&self) -> HashMap<String, f64> {
        let mut d = HashMap::new();
        d.insert("total_variants".into(), self.total_variants as f64);
        d.insert("snps".into(), self.snps as f64);
        d.insert("indels".into(), self.indels as f64);
        d.insert("insertions".into(), self.insertions as f64);
        d.insert("deletions".into(), self.deletions as f64);
        d.insert("pass_filter".into(), self.pass_filter as f64);
        d.insert("ti_tv_ratio".into(), self.ti_tv_ratio);
        d.insert("het_hom_ratio".into(), self.het_hom_ratio);
        d
    }
}

/// VCF file reader.
#[pyclass]
pub struct VcfReader {
    path: String,
}

#[pymethods]
impl VcfReader {
    #[new]
    fn new(path: String) -> PyResult<Self> {
        if !Path::new(&path).exists() {
            return Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
                format!("VCF file not found: {}", path),
            ));
        }
        Ok(VcfReader { path })
    }

    /// Get sample names from VCF header.
    fn samples(&self) -> PyResult<Vec<String>> {
        let reader = open_vcf_reader(&self.path)?;
        
        for line_result in reader.lines() {
            let line = line_result
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
            
            if line.starts_with("#CHROM") {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() > 9 {
                    return Ok(parts[9..].iter().map(|s| s.to_string()).collect());
                }
                break;
            }
        }
        Ok(vec![])
    }

    /// Count total variants.
    fn count(&self) -> PyResult<u64> {
        count_variants(&self.path)
    }

    /// Get VCF statistics.
    fn stats(&self) -> PyResult<VcfStats> {
        calculate_vcf_stats(&self.path)
    }

    /// Read all variants into a list.
    fn read_all(&self) -> PyResult<Vec<Variant>> {
        read_variants(&self.path, None, None, None)
    }

    /// Query variants in a region.
    fn query(&self, chrom: &str, start: u64, end: u64) -> PyResult<Vec<Variant>> {
        read_variants(&self.path, Some(chrom), Some(start), Some(end))
    }

    /// Filter variants by quality threshold.
    fn filter_by_quality(&self, min_qual: f64) -> PyResult<Vec<Variant>> {
        let variants = read_variants(&self.path, None, None, None)?;
        Ok(variants
            .into_iter()
            .filter(|v| v.qual.unwrap_or(0.0) >= min_qual)
            .collect())
    }

    /// Get variants that pass filters.
    fn get_pass_variants(&self) -> PyResult<Vec<Variant>> {
        let variants = read_variants(&self.path, None, None, None)?;
        Ok(variants.into_iter().filter(|v| v.is_pass()).collect())
    }

    /// Get the file path.
    fn get_path(&self) -> &str {
        &self.path
    }
}

/// Open VCF file (handles .gz compression).
fn open_vcf_reader(path: &str) -> PyResult<Box<dyn BufRead>> {
    let file = File::open(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    if path.ends_with(".gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Parse a VCF data line into a Variant.
fn parse_vcf_line(line: &str) -> Option<Variant> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 8 {
        return None;
    }

    let chrom = parts[0].to_string();
    let pos: u64 = parts[1].parse().ok()?;
    let id = if parts[2] == "." {
        None
    } else {
        Some(parts[2].to_string())
    };
    let ref_allele = parts[3].to_string();
    let alt_allele = parts[4].to_string();
    let qual = parts[5].parse().ok();
    let filter = parts[6].to_string();

    // Parse INFO field
    let mut info = HashMap::new();
    for item in parts[7].split(';') {
        if let Some((key, value)) = item.split_once('=') {
            info.insert(key.to_string(), value.to_string());
        } else {
            info.insert(item.to_string(), "true".to_string());
        }
    }

    // Parse genotype if present
    let genotype = if parts.len() > 9 {
        let format_fields: Vec<&str> = parts[8].split(':').collect();
        let sample_fields: Vec<&str> = parts[9].split(':').collect();
        
        format_fields
            .iter()
            .position(|&f| f == "GT")
            .and_then(|idx| sample_fields.get(idx).map(|s| s.to_string()))
    } else {
        None
    };

    Some(Variant {
        chrom,
        pos,
        id,
        ref_allele,
        alt_allele,
        qual,
        filter,
        info,
        genotype,
    })
}

/// Count variants in a VCF file.
#[pyfunction]
pub fn count_variants(path: &str) -> PyResult<u64> {
    let reader = open_vcf_reader(path)?;

    let count = reader
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.starts_with('#'))
        .count() as u64;

    Ok(count)
}

/// Read variants from VCF with optional region filter.
fn read_variants(
    path: &str,
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
) -> PyResult<Vec<Variant>> {
    let reader = open_vcf_reader(path)?;
    let mut variants = Vec::new();

    for line_result in reader.lines() {
        let line = line_result
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        if line.starts_with('#') {
            continue;
        }

        if let Some(variant) = parse_vcf_line(&line) {
            // Apply region filter if specified
            if let Some(target_chrom) = chrom {
                if variant.chrom != target_chrom {
                    continue;
                }
                if let (Some(s), Some(e)) = (start, end) {
                    if variant.pos < s || variant.pos > e {
                        continue;
                    }
                }
            }
            variants.push(variant);
        }
    }

    Ok(variants)
}

/// Calculate comprehensive VCF statistics.
#[pyfunction]
pub fn calculate_vcf_stats(path: &str) -> PyResult<VcfStats> {
    let reader = open_vcf_reader(path)?;
    let mut stats = VcfStats::default();
    let mut total_qual: f64 = 0.0;
    let mut qual_count: u64 = 0;
    let mut het_count: u64 = 0;
    let mut hom_count: u64 = 0;

    for line_result in reader.lines() {
        let line = line_result
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        if line.starts_with('#') {
            continue;
        }

        if let Some(variant) = parse_vcf_line(&line) {
            stats.total_variants += 1;

            // Variant type
            if variant.is_snp() {
                stats.snps += 1;

                // Ti/Tv calculation for SNPs
                let is_transition = matches!(
                    (variant.ref_allele.as_str(), variant.alt_allele.as_str()),
                    ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C")
                );
                if is_transition {
                    stats.transitions += 1;
                } else {
                    stats.transversions += 1;
                }
            } else if variant.is_indel() {
                stats.indels += 1;
                if variant.is_insertion() {
                    stats.insertions += 1;
                } else {
                    stats.deletions += 1;
                }
            }

            // Multi-allelic check
            if variant.alt_allele.contains(',') {
                stats.multi_allelic += 1;
            }

            // Filter status
            if variant.is_pass() {
                stats.pass_filter += 1;
            } else {
                stats.filtered += 1;
            }

            // Quality
            if let Some(q) = variant.qual {
                total_qual += q;
                qual_count += 1;
            }

            // Zygosity
            if variant.is_heterozygous() {
                het_count += 1;
            } else if variant.is_homozygous() {
                hom_count += 1;
            }
        }
    }

    // Calculate ratios
    if stats.transversions > 0 {
        stats.ti_tv_ratio = stats.transitions as f64 / stats.transversions as f64;
    }

    if hom_count > 0 {
        stats.het_hom_ratio = het_count as f64 / hom_count as f64;
    }

    if qual_count > 0 {
        stats.mean_quality = total_qual / qual_count as f64;
    }

    Ok(stats)
}

/// Filter variants by quality and return count.
#[pyfunction]
pub fn count_high_quality_variants(path: &str, min_qual: f64) -> PyResult<u64> {
    let variants = read_variants(path, None, None, None)?;
    let count = variants
        .iter()
        .filter(|v| v.qual.unwrap_or(0.0) >= min_qual)
        .count() as u64;
    Ok(count)
}
