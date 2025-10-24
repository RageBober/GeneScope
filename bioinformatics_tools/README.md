# 🔬 Bioinformatics Tools Integration

## Architecture

```
bioinformatics_tools/
├── base_tool.py          # Abstract base class for all tools
├── config.py             # Configuration management
├── validators.py         # Input validation & security
├── metagenomics/         # Metagenomic analysis tools
│   ├── metagraph.py      # MetaGraph sequence search
│   ├── kraken2.py        # Kraken2 taxonomic classification
│   └── taxonomy.py       # Taxonomy utilities
├── assembly/             # Genome assembly tools
│   ├── megahit.py        # MEGAHIT assembler
│   └── spades.py         # SPAdes assembler
├── classification/       # Classification tools
│   └── gtdb_tk.py        # GTDB-Tk bacterial classifier
├── ml_models/            # Machine learning models
│   ├── xtrimopglm.py     # XtriMoPGLM protein prediction
│   └── enformer.py       # Enformer mutation effects
└── docker/               # Docker configurations
    ├── Dockerfile.kraken2
    ├── Dockerfile.megahit
    └── docker-compose.yml
```

## Quick Start

### Installation

```bash
# Install with bioinformatics tools support
poetry install --extras "biotools io gff"

# Or using pip
pip install -e ".[biotools]"
```

### Basic Usage

```python
from bioinformatics_tools.metagenomics import classify_kraken2

# Classify metagenomic sample
results = classify_kraken2(
    input_path="sample.fastq",
    database_path="/data/kraken2_std",
    confidence=0.1
)

print(f"Classified: {results['results']['classified_percentage']:.1f}%")
```

## Security Model

### Defense-in-Depth Layers

1. **Input Validation**
   - Path traversal prevention
   - File extension whitelisting
   - Size limits
   - Sequence validation

2. **Container Isolation**
   - Docker containerization
   - Non-root users
   - Read-only root filesystem
   - No network access
   - Resource limits

3. **Command Injection Prevention**
   - Argument sanitization
   - No shell execution
   - Subprocess with list arguments

### Example: Secure Tool Execution

```python
from bioinformatics_tools.metagenomics import Kraken2Tool
from bioinformatics_tools.config import Kraken2Config, DockerConfig

# Configure security settings
docker_config = DockerConfig(
    cpu_count=8,
    memory_limit="64g",
    read_only_root=True,
    no_new_privileges=True,
    user="biouser",
    network_mode="none",
    timeout_seconds=3600
)

config = Kraken2Config(
    database_path="/data/kraken2_db",
    docker_config=docker_config
)

tool = Kraken2Tool(config)

# All inputs are validated automatically
try:
    results = tool.classify("sample.fastq")
except ValidationError as exc:
    print(f"Invalid input: {exc}")
```

## Tool Reference

### MetaGraph

**Purpose**: DNA/RNA sequence indexing and search

```python
from bioinformatics_tools.metagenomics import search_metagraph

results = search_metagraph(
    "query.fasta",
    index_dir="/data/metagraph_index",
    discovery_fraction=0.8
)
```

### Kraken2

**Purpose**: Taxonomic classification

```python
from bioinformatics_tools.metagenomics import classify_kraken2, get_top_taxa

results = classify_kraken2(
    "sample.fastq",
    database_path="/data/kraken2_std"
)

top_species = get_top_taxa(results, n=10, rank='S')
```

### MEGAHIT

**Purpose**: Genome assembly

```python
from bioinformatics_tools.assembly import assemble_megahit

results = assemble_megahit(
    "reads_R1.fastq",
    r2_path="reads_R2.fastq",
    preset="meta-sensitive"
)

print(f"N50: {results['results']['n50']} bp")
```

### SPAdes

**Purpose**: Versatile genome assembly

```python
from bioinformatics_tools.assembly import SPAdesTool
from bioinformatics_tools.config import SPAdesConfig

config = SPAdesConfig(mode="meta", careful_mode=True)
tool = SPAdesTool(config)

results = tool.assemble("reads.fastq")
```

### GTDB-Tk

**Purpose**: Bacterial/archaeal classification

```python
from bioinformatics_tools.classification import GTDBTkTool
from bioinformatics_tools.config import GTDBTkConfig

config = GTDBTkConfig(database_path="/data/gtdbtk_db")
tool = GTDBTkTool(config)

results = tool.classify("/genomes")
```

## Configuration

### Environment Variables

```bash
export KRAKEN2_DB_PATH=/data/kraken2_db
export GTDBTK_DATA_PATH=/data/gtdbtk_db
export GENOSCOPE_WORK_DIR=/tmp/genoscope
```

### Config File

Create `~/.config/genoscope/tools.yml`:

```yaml
defaults:
  use_docker: true
  work_dir: "/tmp/genoscope"
  cache_dir: "~/.cache/genoscope"

kraken2:
  database_path: "/data/kraken2_std"
  confidence_threshold: 0.1
  docker:
    memory_limit: "64g"
    cpu_count: 8

megahit:
  min_contig_length: 500
  docker:
    memory_limit: "128g"
    cpu_count: 16
    timeout_seconds: 86400
```

## Docker Setup

### Using Docker Compose

```bash
cd bioinformatics_tools/docker

# Set environment
export KRAKEN2_DB_PATH=/data/kraken2_db
export INPUT_DIR=./input
export OUTPUT_DIR=./output

# Start services
docker-compose up -d

# Run analysis
docker-compose run kraken2 \
  kraken2 --db /db --output /output/results.txt /input/sample.fastq
```

### Building Images

```bash
# Build all images
docker-compose build

# Or build individually
docker build -f Dockerfile.kraken2 -t genoscope/kraken2:2.1.3 .
docker build -f Dockerfile.megahit -t genoscope/megahit:1.2.9 .
```

## Testing

```bash
# Run tests
pytest tests/bioinformatics_tools/

# With coverage
pytest tests/bioinformatics_tools/ --cov=bioinformatics_tools

# Specific test
pytest tests/bioinformatics_tools/test_validators.py::TestValidateFilePath
```

## Performance Benchmarks

| Tool | Input Size | CPU | Memory | Time |
|------|-----------|-----|--------|------|
| Kraken2 | 10M reads | 8 cores | 50GB | ~10 min |
| MEGAHIT | 10M reads | 16 cores | 100GB | ~2 hours |
| SPAdes | 5M reads | 16 cores | 128GB | ~4 hours |
| GTDB-Tk | 100 genomes | 16 cores | 120GB | ~8 hours |

## Troubleshooting

### Docker Permission Denied

```bash
sudo usermod -aG docker $USER
newgrp docker
```

### Database Not Found

```python
# Verify database path
from pathlib import Path
db_path = Path("/data/kraken2_db")
assert db_path.exists(), f"Database not found: {db_path}"
```

### Out of Memory

```python
# Increase memory limit
from bioinformatics_tools.config import DockerConfig

docker_config = DockerConfig(
    memory_limit="256g",
    memory_swap_limit="512g"
)
```

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for development guidelines.

## License

See [LICENSE](../../LICENSE) file.
