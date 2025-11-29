//! BioForge Core - High-performance bioinformatics library
//!
//! This Rust library provides high-performance implementations for:
//! - FASTQ parsing and QC (needletail)
//! - BAM/SAM operations (noodles)
//! - VCF parsing and filtering
//!
//! Exposed to Python via PyO3.

use pyo3::prelude::*;

mod alignment;
mod fastq;
mod utils;
mod variants;

/// BioForge Core Python module.
///
/// Usage from Python:
/// ```python
/// import bioforge_core
///
/// # FASTQ operations
/// stats = bioforge_core.calculate_stats("sample.fastq.gz")
/// print(f"Reads: {stats.total_reads}, GC: {stats.gc_content:.1f}%")
///
/// # BAM operations
/// bam = bioforge_core.BamReader("aligned.bam")
/// aln_stats = bam.stats()
/// print(f"Mapping rate: {aln_stats.mapping_rate:.1f}%")
///
/// # VCF operations
/// vcf = bioforge_core.VcfReader("variants.vcf.gz")
/// vcf_stats = vcf.stats()
/// print(f"Ti/Tv ratio: {vcf_stats.ti_tv_ratio:.2f}")
/// ```
#[pymodule]
fn bioforge_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Version info
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // ─────────────────────────────────────────────────────────────────────────
    // FASTQ module - High-performance FASTQ parsing
    // ─────────────────────────────────────────────────────────────────────────
    m.add_class::<fastq::FastqReader>()?;
    m.add_class::<fastq::FastqRecord>()?;
    m.add_class::<fastq::FastqStats>()?;
    m.add_function(wrap_pyfunction!(fastq::count_reads, m)?)?;
    m.add_function(wrap_pyfunction!(fastq::calculate_stats, m)?)?;

    // ─────────────────────────────────────────────────────────────────────────
    // Alignment module - BAM/SAM operations
    // ─────────────────────────────────────────────────────────────────────────
    m.add_class::<alignment::BamReader>()?;
    m.add_class::<alignment::AlignmentStats>()?;
    m.add_function(wrap_pyfunction!(alignment::get_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::mean_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(alignment::count_mapped_reads, m)?)?;

    // ─────────────────────────────────────────────────────────────────────────
    // Variants module - VCF parsing and statistics
    // ─────────────────────────────────────────────────────────────────────────
    m.add_class::<variants::VcfReader>()?;
    m.add_class::<variants::Variant>()?;
    m.add_class::<variants::VcfStats>()?;
    m.add_function(wrap_pyfunction!(variants::count_variants, m)?)?;
    m.add_function(wrap_pyfunction!(variants::calculate_vcf_stats, m)?)?;
    m.add_function(wrap_pyfunction!(variants::count_high_quality_variants, m)?)?;

    Ok(())
}
