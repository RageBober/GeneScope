# BioForge Core

High-performance Rust library for bioinformatics operations, exposed to Python via PyO3.

## Features

- **FASTQ Processing**: Fast parsing and QC statistics using needletail
- **BAM/SAM Operations**: Alignment statistics and coverage calculation using noodles
- **VCF Parsing**: Variant counting, filtering, and statistics

## Building

### Prerequisites

- Rust toolchain (1.70+)
- Python 3.11+ with development headers
- maturin (`pip install maturin`)

### Development Build

```bash
# Set compatibility flag for Python 3.14+
export PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1  # Linux/Mac
$env:PYO3_USE_ABI3_FORWARD_COMPATIBILITY="1"  # PowerShell

# Build in development mode
cd bioforge_core
maturin develop --release
```

### Release Build

```bash
maturin build --release
pip install target/wheels/bioforge_core-*.whl
```

## Usage

```python
import bioforge_core

# FASTQ statistics
stats = bioforge_core.calculate_stats("sample.fastq.gz")
print(f"Total reads: {stats.total_reads}")
print(f"Mean quality: {stats.mean_quality:.1f}")
print(f"GC content: {stats.gc_content:.1f}%")
print(f"Q30 bases: {stats.q30_bases}")

# Alternative: use FastqReader class
reader = bioforge_core.FastqReader("sample.fastq")
stats = reader.stats()
records = reader.read_all()  # Load all records

# BAM statistics
bam = bioforge_core.BamReader("aligned.bam")
aln_stats = bam.stats()
print(f"Total reads: {aln_stats.total_reads}")
print(f"Mapped: {aln_stats.mapped_reads}")
print(f"Mapping rate: {aln_stats.mapping_rate:.1f}%")
print(f"Mean MAPQ: {aln_stats.mean_mapq:.1f}")

# Coverage for a region
coverage = bam.coverage("chr1", 100000, 100500)

# VCF statistics
vcf = bioforge_core.VcfReader("variants.vcf.gz")
vcf_stats = vcf.stats()
print(f"Total variants: {vcf_stats.total_variants}")
print(f"SNPs: {vcf_stats.snps}")
print(f"Indels: {vcf_stats.indels}")
print(f"Ti/Tv ratio: {vcf_stats.ti_tv_ratio:.2f}")

# Get variants
variants = vcf.read_all()
pass_variants = vcf.get_pass_variants()
region_variants = vcf.query("chr1", 1000000, 2000000)
```

## API Reference

### FASTQ Module

- `FastqReader(path)` - FASTQ file reader
  - `.stats()` -> FastqStats
  - `.read_all()` -> List[FastqRecord]
- `FastqRecord` - Single FASTQ record
  - `.id`, `.sequence`, `.quality`
  - `.mean_quality()`, `.gc_content()`, `.len()`
- `FastqStats` - QC statistics
  - `.total_reads`, `.total_bases`, `.mean_length`
  - `.gc_content`, `.mean_quality`, `.q20_bases`, `.q30_bases`
- `count_reads(path)` -> int
- `calculate_stats(path)` -> FastqStats

### Alignment Module

- `BamReader(path)` - BAM file reader
  - `.stats()` -> AlignmentStats
  - `.coverage(chrom, start, end)` -> List[int]
  - `.reads_per_chromosome()` -> Dict[str, int]
  - `.reference_sequences()` -> List[Tuple[str, int]]
- `AlignmentStats` - Alignment statistics
  - `.total_reads`, `.mapped_reads`, `.unmapped_reads`
  - `.properly_paired`, `.duplicates`
  - `.mapping_rate`, `.mean_mapq`, `.mean_insert_size`
- `get_coverage(bam_path, chrom, start, end)` -> List[int]
- `count_mapped_reads(bam_path)` -> int

### Variants Module

- `VcfReader(path)` - VCF file reader
  - `.stats()` -> VcfStats
  - `.count()` -> int
  - `.samples()` -> List[str]
  - `.read_all()` -> List[Variant]
  - `.query(chrom, start, end)` -> List[Variant]
  - `.get_pass_variants()` -> List[Variant]
  - `.filter_by_quality(min_qual)` -> List[Variant]
- `Variant` - Single variant
  - `.chrom`, `.pos`, `.ref_allele`, `.alt_allele`
  - `.qual`, `.filter`, `.genotype`
  - `.is_snp()`, `.is_indel()`, `.variant_type()`
  - `.is_pass()`, `.is_homozygous()`, `.is_heterozygous()`
- `VcfStats` - VCF statistics
  - `.total_variants`, `.snps`, `.indels`
  - `.ti_tv_ratio`, `.het_hom_ratio`, `.mean_quality`
- `count_variants(path)` -> int
- `calculate_vcf_stats(path)` -> VcfStats

## Performance

Benchmarks on typical WGS data:

| Operation | File Size | Time |
|-----------|-----------|------|
| FASTQ stats | 10GB | ~30s |
| BAM stats | 50GB | ~2min |
| VCF count | 1M variants | ~5s |

## License

MIT
