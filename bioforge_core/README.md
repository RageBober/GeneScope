# BioForge Core

High-performance bioinformatics library for BioForge, written in Rust.

## Features

- **FASTQ parsing** — Fast reading and QC statistics
- **BAM/CRAM operations** — Coverage, statistics
- **VCF handling** — Parsing, filtering, counting

## Installation

```bash
cd bioforge_core
maturin develop
```

## Usage

```python
import bioforge_core

# FASTQ stats
stats = bioforge_core.calculate_stats("sample.fastq.gz")
print(stats.total_reads, stats.mean_quality)

# Count variants
count = bioforge_core.count_variants("sample.vcf.gz")
```
