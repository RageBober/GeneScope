//! VCF parsing and filtering module.
//!
//! High-performance VCF operations using noodles.

use pyo3::prelude::*;
use std::collections::HashMap;

/// Variant representation.
#[pyclass]
#[derive(Clone)]
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
}

#[pymethods]
impl Variant {
    #[new]
    fn new(
        chrom: String,
        pos: u64,
        ref_allele: String,
        alt_allele: String,
    ) -> Self {
        Variant {
            chrom,
            pos,
            id: None,
            ref_allele,
            alt_allele,
            qual: None,
            filter: "PASS".to_string(),
            info: HashMap::new(),
        }
    }
    
    fn __repr__(&self) -> String {
        format!(
            "Variant({}:{} {}>{}, filter={})",
            self.chrom, self.pos, self.ref_allele, self.alt_allele, self.filter
        )
    }
    
    /// Check if variant is a SNP.
    fn is_snp(&self) -> bool {
        self.ref_allele.len() == 1 && self.alt_allele.len() == 1
    }
    
    /// Check if variant is an indel.
    fn is_indel(&self) -> bool {
        self.ref_allele.len() != self.alt_allele.len()
    }
    
    /// Get variant type.
    fn variant_type(&self) -> String {
        if self.is_snp() {
            "SNV".to_string()
        } else if self.is_indel() {
            if self.ref_allele.len() > self.alt_allele.len() {
                "DEL".to_string()
            } else {
                "INS".to_string()
            }
        } else {
            "MNV".to_string()
        }
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
        Ok(VcfReader { path })
    }
    
    /// Get samples in VCF.
    fn samples(&self) -> PyResult<Vec<String>> {
        // TODO: Implement using noodles-vcf
        Ok(vec![])
    }
    
    /// Count variants.
    fn count(&self) -> PyResult<u64> {
        count_variants(&self.path)
    }
    
    /// Get variants in a region.
    fn query(&self, _chrom: &str, _start: u64, _end: u64) -> PyResult<Vec<Variant>> {
        // TODO: Implement region query with tabix
        Ok(vec![])
    }
    
    /// Get the file path.
    fn get_path(&self) -> &str {
        &self.path
    }
}

/// Count variants in a VCF file.
#[pyfunction]
pub fn count_variants(path: &str) -> PyResult<u64> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use flate2::read::GzDecoder;
    
    let file = File::open(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    
    let count = reader
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.starts_with('#'))
        .count() as u64;
    
    Ok(count)
}
