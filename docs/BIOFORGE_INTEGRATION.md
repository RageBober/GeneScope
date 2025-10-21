# 🔬 BioForge Integration Guide

## Overview

This guide explains how GeneScope integrates external bioinformatics tools into unified **BioForge workflows**. Each tool serves a specific purpose in the genomic analysis pipeline, from data discovery to mutation effect prediction.

---

## 🧬 Tool Ecosystem & Roles

### Summary Table

| Инструмент | Описание | Роль в BioForge |
|------------|----------|----------------|
| **MetaGraph** | Индексирует и ищет последовательности ДНК/РНК по всему миру | Источник данных и связка с GeneScope |
| **Kraken2 / Centrifuge** | Таксономическая классификация метагеномных данных (находит, какой организм присутствует) | Автоматическое определение происхождения найденных генов |
| **MEGAHIT / SPAdes** | Быстрая сборка геномов из коротких ридов | Встроенный сборщик для входных FASTQ-файлов |
| **GTDB-Tk** | Классификация бактерий по геномной базе GTDB | Нормализация и унификация таксономических меток |
| **XtriMoPGLM** | Модель Transformer для предсказания белковых свойств | Аналитика аминокислотных последовательностей перед симуляцией |
| **DeepSEA / Enformer** | Предсказывают влияние мутаций ДНК на экспрессию генов | Использовать в EvoSim для прогноза фенотипических эффектов |

---

## 📊 BioForge Pipelines

### 1. MetagenomicsPipeline

**Purpose**: End-to-end metagenomic sample analysis

**Workflow**:
```
FASTQ reads
    ↓
[Kraken2] → Taxonomic Classification
    ↓
[MEGAHIT/SPAdes] → Genome Assembly
    ↓
[GTDB-Tk] → Taxonomy Normalization
    ↓
Classified & Assembled Genomes
```

**Code Example**:
```python
from bioinformatics_tools import MetagenomicsPipeline

# Initialize pipeline
pipeline = MetagenomicsPipeline(
    kraken2_db="/data/kraken2_standard",
    gtdbtk_db="/data/gtdbtk_r214",
    assembler="megahit"
)

# Run complete analysis
result = pipeline.run(
    fastq_r1="gut_sample_R1.fastq.gz",
    fastq_r2="gut_sample_R2.fastq.gz",
    output_dir="results/metagenomics"
)

# Access results
print(f"Classified: {result.results['kraken2']['classified_percentage']:.1f}%")
print(f"Assembly N50: {result.results['assembly']['n50']} bp")
print(f"Top species: {result.results['kraken2_top_species'][0]['name']}")
```

**Use Cases**:
- 🦠 **Microbiome Studies**: Analyze human gut, soil, ocean samples
- 🏥 **Clinical Diagnostics**: Identify pathogens in patient samples
- 🌍 **Environmental Genomics**: Study microbial communities in ecosystems

---

### 2. ProteinAnalysisPipeline

**Purpose**: Predict protein properties using AI/ML models

**Workflow**:
```
Protein FASTA
    ↓
[XtriMoPGLM Transformer] → Property Prediction
    ↓
Function, Structure, Stability predictions
```

**Code Example**:
```python
from bioinformatics_tools import ProteinAnalysisPipeline

# Initialize with GPU
pipeline = ProteinAnalysisPipeline(
    device="cuda",
    batch_size=64
)

# Analyze proteins
result = pipeline.run(
    protein_fasta="novel_proteins.fasta",
    output_dir="results/proteins"
)

# Predictions available as JSON
predictions = result.results['predictions']
```

**Use Cases**:
- 🧪 **Protein Engineering**: Optimize enzyme activity before lab synthesis
- 💊 **Drug Discovery**: Predict drug target interactions
- 🔬 **Functional Annotation**: Characterize unknown proteins from metagenomes

**Integration with EvoSim**:
```python
# Use predictions in evolutionary simulations
from evosim import ProteinEvolution  # hypothetical

evolution = ProteinEvolution(
    initial_sequence="MKTAYIAKQR...",
    fitness_predictor=lambda seq: xtrimopglm.predict(seq)['stability']
)

evolved_proteins = evolution.run(generations=1000)
```

---

### 3. MutationEffectPipeline

**Purpose**: Predict phenotypic effects of DNA mutations

**Workflow**:
```
DNA sequences (with variants)
    ↓
[Enformer/DeepSEA] → Expression Prediction
    ↓
Mutation effect scores
```

**Code Example**:
```python
from bioinformatics_tools import MutationEffectPipeline

# Initialize Enformer model
pipeline = MutationEffectPipeline(
    model_name="enformer",
    device="cuda"
)

# Predict effects
result = pipeline.run(
    sequence_fasta="candidate_mutations.fasta",
    output_dir="results/mutations"
)

# Effect predictions for each mutation
effects = result.results['predictions']
```

**Use Cases**:
- 🧬 **Variant Interpretation**: Understand disease-causing mutations
- 🌱 **Crop Improvement**: Predict effects of breeding targets
- 🦠 **Pathogen Evolution**: Forecast adaptive mutations

**EvoSim Integration**:
```python
# Use mutation effects to guide evolutionary simulations
from evosim import AdaptiveEvolution

# Define fitness landscape based on Enformer predictions
def fitness_function(genome):
    mutation_effects = enformer.predict(genome)
    return mutation_effects['expression_score']

simulation = AdaptiveEvolution(
    population_size=1000,
    mutation_rate=1e-6,
    fitness_func=fitness_function
)

final_population = simulation.run(generations=10000)
```

---

### 4. MetaGraphDataSource

**Purpose**: Global sequence discovery and data integration

**Workflow**:
```
Query sequences
    ↓
[MetaGraph Index Search]
    ↓
Matches in SRA/public databases
    ↓
Download & integrate into GeneScope
```

**Code Example**:
```python
from bioinformatics_tools import MetaGraphDataSource

# Connect to MetaGraph index
data_source = MetaGraphDataSource(
    index_dir="/data/metagraph_sra_index"
)

# Search for homologous sequences
result = data_source.search_sequences(
    query_fasta="gene_of_interest.fasta",
    discovery_fraction=0.8
)

print(f"Found {result.results['num_matches']} matches in public databases")

# Matches can be downloaded and fed into other pipelines
for match in result.results['matches']:
    print(f"  - {match['matched_labels']}")
```

**Use Cases**:
- 📚 **Data Mining**: Find related sequences in massive public datasets
- 🔗 **Cross-Study Integration**: Connect findings across multiple studies
- 🌐 **Global Pathogen Surveillance**: Track variants worldwide

**GeneScope Integration**:
```python
# Use MetaGraph to enrich GeneScope analysis
from data_analysis import load_data

# Load initial gene set
genes = load_data("initial_genes.csv", "csv")

# Search for homologs in MetaGraph
metagraph_results = data_source.search_sequences("initial_genes.fasta")

# Integrate found sequences
enriched_dataset = merge_metagraph_results(genes, metagraph_results)

# Continue with GeneScope analysis
from data_analysis import filter_outliers, extract_pca
filtered = filter_outliers(enriched_dataset, method="iqr")
pca_result = extract_pca(filtered)
```

---

## 🔄 Multi-Tool Workflows

### Complete Metagenomics + Protein Characterization

```python
from bioinformatics_tools import (
    MetagenomicsPipeline,
    ProteinAnalysisPipeline,
    MutationEffectPipeline
)

# Stage 1: Metagenomics
mg_pipeline = MetagenomicsPipeline(
    kraken2_db="/data/kraken2_standard",
    assembler="megahit"
)

mg_result = mg_pipeline.run(
    fastq_r1="sample_R1.fastq.gz",
    fastq_r2="sample_R2.fastq.gz",
    output_dir="results/metagenomics"
)

# Stage 2: Extract proteins from contigs
# (Use tool like Prodigal to predict genes)
# prodigal -i contigs.fasta -a proteins.fasta

# Stage 3: Protein analysis
protein_pipeline = ProteinAnalysisPipeline(device="cuda")

protein_result = protein_pipeline.run(
    protein_fasta="results/metagenomics/proteins.fasta",
    output_dir="results/proteins"
)

# Stage 4: Mutation effect prediction
mutation_pipeline = MutationEffectPipeline(model_name="enformer")

mutation_result = mutation_pipeline.run(
    sequence_fasta="results/metagenomics/variant_regions.fasta",
    output_dir="results/mutations"
)

# Combined results
print(f"Organisms found: {len(mg_result.results['kraken2_top_species'])}")
print(f"Proteins analyzed: {len(protein_result.results['predictions'])}")
print(f"Mutations scored: {len(mutation_result.results['predictions'])}")
```

---

## 🔗 Integration with Existing GeneScope Modules

### Data Flow Diagram

```
External Data Sources
    ↓
[MetaGraph] → Search & Download
    ↓
GeneScope Data Ingestion (data_analysis/data_ingestion.py)
    ↓
[load_data] → CSV/FASTA/VCF/BAM
    ↓
Data Cleaning (data_analysis/data_cleaning.py)
    ↓
BioForge Pipelines
    ├→ [MetagenomicsPipeline] → Taxonomic insights
    ├→ [ProteinAnalysisPipeline] → Functional predictions
    └→ [MutationEffectPipeline] → Phenotype forecasts
    ↓
GeneScope Analysis (data_analysis/analysis_core.py)
    ├→ Statistical analysis
    ├→ PCA / dimensionality reduction
    └→ Feature selection
    ↓
Visualization (data_analysis/visualization.py)
    └→ Plots, heatmaps, PCA graphs
    ↓
Export Results
```

### Example: Integrated Analysis

```python
# Step 1: Load data from various sources
from data_analysis import load_data
from bioinformatics_tools import MetaGraphDataSource

# Load existing dataset
genes = load_data("genes.csv", "csv")

# Enrich with MetaGraph
metagraph = MetaGraphDataSource(index_dir="/data/metagraph_index")
homologs = metagraph.search_sequences("genes.fasta")

# Step 2: Run metagenomics pipeline
from bioinformatics_tools import MetagenomicsPipeline

mg_pipeline = MetagenomicsPipeline(kraken2_db="/data/kraken2_db")
mg_result = mg_pipeline.run("sample.fastq")

# Step 3: Clean and filter data
from data_analysis import remove_duplicates, handle_missing_values, filter_outliers

cleaned = remove_duplicates(genes)
cleaned = handle_missing_values(cleaned, method="mean")
filtered = filter_outliers(cleaned, "expression_level", method="iqr")

# Step 4: Protein prediction
from bioinformatics_tools import ProteinAnalysisPipeline

protein_pipeline = ProteinAnalysisPipeline()
protein_predictions = protein_pipeline.run("proteins.fasta")

# Step 5: Statistical analysis
from data_analysis import extract_pca, calculate_correlation

pca_result = extract_pca(filtered, n_components=3)
correlations = calculate_correlation(filtered)

# Step 6: Visualization
from data_analysis import plot_pca, plot_heatmap

plot_pca(pca_result, target_column="organism")
plot_heatmap(correlations)

print("✅ Integrated analysis complete!")
```

---

## 🔒 Security & Performance

### Security Features (All Pipelines)

1. **Input Validation**
   - Path traversal prevention
   - File type whitelisting
   - Sequence validation (IUPAC codes)

2. **Container Isolation**
   - Docker containerization
   - Non-root execution
   - Read-only filesystems
   - No network access

3. **Resource Limiting**
   - CPU quotas
   - Memory limits
   - Execution timeouts

### Performance Benchmarks

| Pipeline | Typical Input | Processing Time | Resource Usage |
|----------|--------------|----------------|---------------|
| **MetagenomicsPipeline** | 10M reads (paired) | 2-3 hours | 8 cores, 64GB RAM |
| **ProteinAnalysisPipeline** | 1000 proteins | 1-2 minutes | 1 GPU, 16GB VRAM |
| **MutationEffectPipeline** | 100 variants | 5-10 minutes | 1 GPU, 32GB VRAM |
| **MetaGraphDataSource** | 100 queries | <1 minute | 4 cores, 8GB RAM |

---

## 📚 References & Citations

### Tools Documentation

- **MetaGraph**: [github.com/ratschlab/metagraph](https://github.com/ratschlab/metagraph)
- **Kraken2**: [github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2)
- **MEGAHIT**: [github.com/voutcn/megahit](https://github.com/voutcn/megahit)
- **SPAdes**: [github.com/ablab/spades](https://github.com/ablab/spades)
- **GTDB-Tk**: [github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)
- **Enformer**: [deepmind.com/research/enformer](https://www.deepmind.com/research)

### Key Publications

1. **Kraken2**: Wood et al. (2019). "Improved metagenomic analysis with Kraken 2." *Genome Biology* 20:257.

2. **MEGAHIT**: Li et al. (2015). "MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly." *Bioinformatics* 31(10):1674-1676.

3. **GTDB**: Parks et al. (2022). "GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy." *Nucleic Acids Research* 50(D1):D785-D794.

4. **Enformer**: Avsec et al. (2021). "Effective gene expression prediction from sequence by integrating long-range interactions." *Nature Methods* 18:1196–1203.

---

## 🚀 Quick Start

### Installation

```bash
# Install GeneScope with BioForge tools
poetry install --extras "biotools io gff"

# Or with pip
pip install -e ".[biotools,io,gff]"
```

### Download Databases

```bash
# Kraken2 standard database (~50GB)
kraken2-build --standard --db /data/kraken2_standard

# GTDB-Tk database (~80GB)
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xzf gtdbtk_data.tar.gz -C /data/gtdbtk_db
```

### Run Example Workflow

```bash
# Navigate to examples
cd examples

# Run BioForge workflows
python bioforge_workflows.py
```

---

## 💡 Best Practices

1. **Start with Dry Run Mode**
   ```python
   config = MetagenomicsPipeline(kraken2_db="/data/kraken2", dry_run=True)
   # This shows commands without execution
   ```

2. **Monitor Resource Usage**
   ```bash
   docker stats  # Monitor container resources
   htop  # Monitor system resources
   ```

3. **Save Intermediate Results**
   ```python
   result.to_json("results/pipeline_result.json")
   ```

4. **Use Appropriate Assembler**
   - **MEGAHIT**: Fast, low memory (metagenomes, drafts)
   - **SPAdes**: Accurate, high memory (isolates, single-cell)

5. **Optimize Batch Sizes**
   ```python
   # For GPU workloads, tune based on VRAM
   protein_pipeline = ProteinAnalysisPipeline(
       batch_size=128 if gpu_memory > 24 else 64
   )
   ```

---

## 🐛 Troubleshooting

### Common Issues

**Issue**: Out of memory during MEGAHIT assembly

**Solution**:
```python
config = MEGAHITConfig(
    docker_config=DockerConfig(
        memory_limit="256g",  # Increase limit
        memory_swap_limit="512g"
    )
)
```

**Issue**: Kraken2 database not found

**Solution**:
```python
# Verify database path
from pathlib import Path
assert Path("/data/kraken2_db").exists()

# Or set environment variable
import os
os.environ['KRAKEN2_DB_PATH'] = "/data/kraken2_db"
```

**Issue**: GPU out of memory for ML models

**Solution**:
```python
# Reduce batch size or use CPU
config = XtriMoPGLMConfig(
    batch_size=16,  # Reduce from 64
    device="cpu"  # Fallback to CPU
)
```

---

## 🎯 Next Steps

1. ✅ Review [examples/bioforge_workflows.py](../examples/bioforge_workflows.py)
2. ✅ Download required databases
3. ✅ Test on sample data
4. ✅ Integrate with existing GeneScope analysis
5. ✅ Deploy in production with Docker Compose

For more details, see:
- [BIOINFORMATICS_TOOLS.md](./BIOINFORMATICS_TOOLS.md) - Tool-specific documentation
- [bioinformatics_tools/README.md](../bioinformatics_tools/README.md) - Module architecture
