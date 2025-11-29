//! BioForge Core - High-performance bioinformatics library
//!
//! This Rust library provides high-performance implementations for:
//! - FASTQ parsing and QC
//! - BAM/CRAM operations
//! - VCF parsing and filtering
//!
//! Exposed to Python via PyO3.

use pyo3::prelude::*;

mod fastq;
mod alignment;
mod variants;
mod utils;

/// BioForge Core Python module.
#[pymodule]
fn bioforge_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Version info
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    
    // FASTQ module
    m.add_class::<fastq::FastqReader>()?;
    m.add_class::<fastq::FastqRecord>()?;
    m.add_class::<fastq::FastqStats>()?;
    m.add_function(wrap_pyfunction!(fastq::count_reads, m)?)?;
    m.add_function(wrap_pyfunction!(fastq::calculate_stats, m)?)?;
    
    // Alignment module
    m.add_class::<alignment::BamReader>()?;
    m.add_class::<alignment::AlignmentStats>()?;
    m.add_function(wrap_pyfunction!(alignment::get_coverage, m)?)?;
    
    // Variants module
    m.add_class::<variants::VcfReader>()?;
    m.add_class::<variants::Variant>()?;
    m.add_function(wrap_pyfunction!(variants::count_variants, m)?)?;
    
    Ok(())
}
