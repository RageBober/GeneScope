# 🧬 Bioinformatics Tools Integration Guide

## Overview

GeneScope integrates several external bioinformatics tools through a secure, containerized architecture. All tools run in isolated Docker containers with resource limits and security controls.

## Integrated Tools

### 1. MetaGraph - DNA/RNA Sequence Indexing

**Purpose**: Index and search large collections of DNA/RNA sequences

**Use Cases**:
- Search for sequences across multiple datasets
- Find similar sequences in large databases
- Metagenomic analysis

**Example Usage**:
```python
from bioinformatics_tools.metagenomics import search_metagraph

# Search for sequences
results = search_metagraph(
    query_path="queries.fasta",
    index_dir="/data/metagraph_index",
    discovery_fraction=0.8
)

print(f"Found {results['results']['num_matches']} matches")
```

### 2. Kraken2 - Taxonomic Classification

**Purpose**: Taxonomic classification of metagenomic samples using exact k-mer matches

**Use Cases**:
- Taxonomic classification of metagenomic samples
- Contamination detection in genomic data
- Microbiome analysis

**Example Usage**:
```python
from bioinformatics_tools.metagenomics import classify_kraken2, get_top_taxa

# Classify sequences
results = classify_kraken2(
    input_path="sample.fastq",
    database_path="/data/kraken2_std",
    confidence=0.1
)

print(f"Classified: {results['results']['classified_percentage']:.1f}%")

# Get top species
top_species = get_top_taxa(results, n=10, rank='S')
for taxon in top_species:
    print(f"{taxon['name']}: {taxon['percentage']:.2f}%")
```

### 3. MEGAHIT - Genome Assembly

**Purpose**: Ultra-fast and memory-efficient NGS assembler

**Use Cases**:
- Metagenomic assembly from short reads
- Low-memory genome assembly
- Fast assembly of complex microbial communities

**Example Usage**:
```python
from bioinformatics_tools.assembly import assemble_megahit

# Assemble paired-end reads
results = assemble_megahit(
    input_path="reads_R1.fastq",
    r2_path="reads_R2.fastq",
    preset="meta-sensitive"
)

print(f"Assembled {results['results']['num_contigs']} contigs")
print(f"N50: {results['results']['n50']} bp")
```

### 4. SPAdes - Versatile Genome Assembler

**Purpose**: Genome assembly for isolates, single-cell, metagenomes, plasmids, RNA-Seq

**Example Usage**:
```python
from bioinformatics_tools.assembly import SPAdesTool
from bioinformatics_tools.config import SPAdesConfig

config = SPAdesConfig(mode="meta", careful_mode=True)
tool = SPAdesTool(config)

results = tool.assemble(
    "reads.fastq",
    output_dir="/results"
)
```

### 5. GTDB-Tk - Bacterial Classification

**Purpose**: Classify bacterial and archaeal genomes using GTDB taxonomy

**Example Usage**:
```python
from bioinformatics_tools.classification import GTDBTkTool
from bioinformatics_tools.config import GTDBTkConfig

config = GTDBTkConfig(
    database_path="/data/gtdbtk_db",
    num_threads=16
)

tool = GTDBTkTool(config)
results = tool.classify(
    genome_dir="/genomes",
    output_dir="/results"
)
```

### 6. XtriMoPGLM - Protein Prediction

**Purpose**: Transformer-based model for protein property prediction

**Example Usage**:
```python
from bioinformatics_tools.ml_models import XtriMoPGLMTool

tool = XtriMoPGLMTool()
results = tool.predict(
    "proteins.fasta",
    output_dir="/results"
)
```

### 7. Enformer - Mutation Effect Prediction

**Purpose**: Predict gene expression and mutation effects from DNA sequence

**Example Usage**:
```python
from bioinformatics_tools.ml_models import EnformerTool

tool = EnformerTool()
results = tool.predict(
    "sequences.fasta",
    output_dir="/results"
)
```

## Security Features

### 1. Input Validation
- Path traversal prevention
- File type whitelisting
- Sequence validation
- Size limits

### 2. Docker Isolation
- All tools run in containers
- Non-root users
- Read-only root filesystems
- No network access by default
- CPU and memory limits

### 3. Resource Management
- Configurable CPU limits
- Memory quotas
- Execution timeouts
- Disk space management

## Configuration

### Global Configuration

Create a configuration file `~/.config/genoscope/config.yml`:

```yaml
bioinformatics_tools:
  # Docker settings
  use_docker: true
  docker:
    cpu_count: 8
    memory_limit: "64g"
    timeout_seconds: 3600

  # Database paths
  kraken2_db: "/data/kraken2_db"
  gtdbtk_db: "/data/gtdbtk_db"

  # Working directories
  work_dir: "/tmp/genoscope"
  cache_dir: "~/.cache/genoscope"
```

### Tool-Specific Configuration

```python
from bioinformatics_tools.config import Kraken2Config, DockerConfig

# Custom Docker configuration
docker_cfg = DockerConfig(
    cpu_count=16,
    memory_limit="128g",
    timeout_seconds=7200,
    user="biouser"
)

# Tool configuration
config = Kraken2Config(
    database_path="/data/kraken2_custom",
    confidence_threshold=0.15,
    docker_config=docker_cfg
)
```

## Docker Deployment

### Using Docker Compose

```bash
# Set environment variables
export KRAKEN2_DB_PATH=/data/kraken2_db
export GTDBTK_DB_PATH=/data/gtdbtk_db
export INPUT_DIR=./input
export OUTPUT_DIR=./output

# Start services
cd bioinformatics_tools/docker
docker-compose up -d

# Run Kraken2
docker-compose run kraken2 \
  kraken2 --db /db --output /output/results.txt /input/sample.fastq
```

### Building Custom Images

```bash
# Build Kraken2 image
docker build -f Dockerfile.kraken2 -t genoscope/kraken2:2.1.3 .

# Build MEGAHIT image
docker build -f Dockerfile.megahit -t genoscope/megahit:1.2.9 .
```

## Best Practices

### 1. Resource Planning

```python
# For large metagenomes
config = MEGAHITConfig(
    docker_config=DockerConfig(
        cpu_count=32,
        memory_limit="256g",
        timeout_seconds=86400  # 24 hours
    )
)
```

### 2. Error Handling

```python
from bioinformatics_tools.base_tool import ToolExecutionError

try:
    results = classify_kraken2("sample.fastq", "/data/kraken2_db")
except ToolExecutionError as exc:
    logger.error(f"Classification failed: {exc}")
    # Handle error
```

### 3. Dry Run Mode

```python
# Test command without execution
config = Kraken2Config(dry_run=True)
tool = Kraken2Tool(config)

results = tool.classify("sample.fastq")
print(results['command'])  # View generated command
```

## Troubleshooting

### Common Issues

**Issue**: Docker permission denied
```bash
# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker
```

**Issue**: Database not found
```python
# Check database path
config = Kraken2Config(database_path="/correct/path")
```

**Issue**: Out of memory
```python
# Increase memory limit
config.docker_config.memory_limit = "128g"
```

## Performance Tips

1. **Use SSD storage** for databases and temporary files
2. **Allocate sufficient RAM** (Kraken2 needs ~50GB for standard DB)
3. **Use tmpfs** for temporary files in Docker containers
4. **Enable parallel processing** where supported
5. **Cache results** for repeated analyses

## References

- [MetaGraph](https://github.com/ratschlab/metagraph)
- [Kraken2](https://github.com/DerrickWood/kraken2)
- [MEGAHIT](https://github.com/voutcn/megahit)
- [SPAdes](https://github.com/ablab/spades)
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)
- [Enformer](https://github.com/deepmind/deepmind-research/tree/master/enformer)
