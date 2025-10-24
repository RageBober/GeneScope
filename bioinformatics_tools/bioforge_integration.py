"""
🔬 BioForge Integration Layer

Orchestrates bioinformatics tools into unified workflows for GeneScope/BioForge.

This module provides high-level pipelines that combine multiple tools:
- Metagenomics: Kraken2 → MEGAHIT/SPAdes → GTDB-Tk
- Protein Analysis: XtriMoPGLM for protein property prediction
- Mutation Effects: Enformer for phenotype prediction
- Data Sources: MetaGraph for sequence discovery
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from .assembly import MEGAHITTool, SPAdesTool
from .classification import GTDBTkTool
from .config import (
    GTDBTkConfig,
    Kraken2Config,
    MEGAHITConfig,
    MetaGraphConfig,
    SPAdesConfig,
    XtriMoPGLMConfig,
    EnformerConfig,
)
from .metagenomics import Kraken2Tool, MetaGraphTool, TaxonomyAnalyzer
from .ml_models import XtriMoPGLMTool, EnformerTool

logger = logging.getLogger(__name__)


@dataclass
class PipelineResult:
    """Result from a BioForge pipeline execution"""

    pipeline_name: str
    status: str  # 'success', 'partial', 'failed'
    steps_completed: list[str] = field(default_factory=list)
    steps_failed: list[str] = field(default_factory=list)
    results: dict[str, Any] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)

    def to_json(self, output_path: Path) -> None:
        """Save pipeline result to JSON file"""
        with open(output_path, 'w') as f:
            json.dump({
                'pipeline_name': self.pipeline_name,
                'status': self.status,
                'steps_completed': self.steps_completed,
                'steps_failed': self.steps_failed,
                'results': self.results,
                'errors': self.errors,
            }, f, indent=2)


# ═══════════════════════════════════════════════════════════════
#  Metagenomics Pipeline
# ═══════════════════════════════════════════════════════════════

class MetagenomicsPipeline:
    """
    End-to-end metagenomics analysis pipeline.

    Workflow:
    1. Taxonomic classification with Kraken2 (identify organisms)
    2. Genome assembly with MEGAHIT or SPAdes
    3. Bacterial classification with GTDB-Tk (normalize taxonomy)

    Use Case:
    - Analyze metagenomic samples (microbiome, environmental)
    - Identify organisms and assemble their genomes
    - Standardize taxonomy using GTDB

    Examples:
        >>> pipeline = MetagenomicsPipeline(
        ...     kraken2_db="/data/kraken2_std",
        ...     gtdbtk_db="/data/gtdbtk_db"
        ... )
        >>> result = pipeline.run(
        ...     fastq_r1="sample_R1.fastq",
        ...     fastq_r2="sample_R2.fastq",
        ...     output_dir="/results"
        ... )
        >>> print(f"Status: {result.status}")
        >>> print(f"Top taxa: {result.results['kraken2_top_taxa']}")
        >>> print(f"Assembly N50: {result.results['assembly_n50']}")
    """

    def __init__(
        self,
        kraken2_db: str | Path,
        gtdbtk_db: Optional[str | Path] = None,
        assembler: str = "megahit",  # 'megahit' or 'spades'
        use_docker: bool = True,
    ):
        """
        Initialize metagenomics pipeline.

        Args:
            kraken2_db: Path to Kraken2 database
            gtdbtk_db: Path to GTDB-Tk database (optional)
            assembler: Which assembler to use ('megahit' or 'spades')
            use_docker: Use Docker containerization
        """
        self.kraken2_db = Path(kraken2_db)
        self.gtdbtk_db = Path(gtdbtk_db) if gtdbtk_db else None
        self.assembler = assembler.lower()
        self.use_docker = use_docker

        # Initialize tools
        self.kraken2 = Kraken2Tool(
            Kraken2Config(
                database_path=self.kraken2_db,
                use_docker=use_docker
            )
        )

        if self.assembler == "megahit":
            self.assembler_tool = MEGAHITTool(
                MEGAHITConfig(use_docker=use_docker)
            )
        elif self.assembler == "spades":
            self.assembler_tool = SPAdesTool(
                SPAdesConfig(mode="meta", use_docker=use_docker)
            )
        else:
            raise ValueError(f"Unknown assembler: {assembler}")

        if self.gtdbtk_db:
            self.gtdbtk = GTDBTkTool(
                GTDBTkConfig(
                    database_path=self.gtdbtk_db,
                    use_docker=use_docker
                )
            )
        else:
            self.gtdbtk = None

        self.taxonomy_analyzer = TaxonomyAnalyzer()

    def run(
        self,
        fastq_r1: str | Path,
        fastq_r2: Optional[str | Path] = None,
        output_dir: Optional[str | Path] = None,
        skip_classification: bool = False,
        skip_assembly: bool = False,
        skip_gtdbtk: bool = False,
    ) -> PipelineResult:
        """
        Run complete metagenomics pipeline.

        Args:
            fastq_r1: Forward reads (or single-end)
            fastq_r2: Reverse reads (for paired-end)
            output_dir: Output directory
            skip_classification: Skip Kraken2 classification
            skip_assembly: Skip genome assembly
            skip_gtdbtk: Skip GTDB-Tk classification

        Returns:
            PipelineResult with all step results
        """
        result = PipelineResult(
            pipeline_name="MetagenomicsPipeline",
            status="running"
        )

        logger.info("Starting metagenomics pipeline")

        # Step 1: Taxonomic Classification with Kraken2
        if not skip_classification:
            try:
                logger.info("Step 1/3: Taxonomic classification (Kraken2)")
                kraken2_result = self.kraken2.classify(
                    fastq_r1,
                    output_dir=output_dir,
                    paired=fastq_r2
                )

                result.steps_completed.append("kraken2_classification")
                result.results["kraken2"] = kraken2_result["results"]

                # Extract top taxa
                top_taxa = self.taxonomy_analyzer.get_top_taxa(
                    kraken2_result["results"].get("taxa_found", []),
                    n=10,
                    rank="S"
                )
                result.results["kraken2_top_species"] = top_taxa

                logger.info(
                    f"Classification complete: "
                    f"{kraken2_result['results']['classified_percentage']:.1f}% classified"
                )

            except Exception as exc:
                logger.error(f"Kraken2 classification failed: {exc}")
                result.steps_failed.append("kraken2_classification")
                result.errors["kraken2"] = str(exc)

        # Step 2: Genome Assembly
        if not skip_assembly:
            try:
                logger.info(f"Step 2/3: Genome assembly ({self.assembler})")
                assembly_result = self.assembler_tool.assemble(
                    fastq_r1,
                    output_dir=output_dir,
                    r2=fastq_r2
                )

                result.steps_completed.append(f"{self.assembler}_assembly")
                result.results["assembly"] = assembly_result["results"]

                logger.info(
                    f"Assembly complete: "
                    f"{assembly_result['results']['num_contigs']} contigs, "
                    f"N50={assembly_result['results']['n50']} bp"
                )

            except Exception as exc:
                logger.error(f"Assembly failed: {exc}")
                result.steps_failed.append(f"{self.assembler}_assembly")
                result.errors["assembly"] = str(exc)

        # Step 3: GTDB-Tk Classification (if contigs available)
        if not skip_gtdbtk and self.gtdbtk and "assembly" in result.results:
            try:
                logger.info("Step 3/3: Bacterial classification (GTDB-Tk)")

                # GTDB-Tk requires genomes in bins (separate directories)
                # For now, we'll skip if no binning was done
                logger.warning(
                    "GTDB-Tk requires genome bins. "
                    "Run binning tool (e.g., MetaBAT) before GTDB-Tk."
                )
                result.steps_completed.append("gtdbtk_skipped")

            except Exception as exc:
                logger.error(f"GTDB-Tk classification failed: {exc}")
                result.steps_failed.append("gtdbtk_classification")
                result.errors["gtdbtk"] = str(exc)

        # Determine final status
        if result.steps_failed:
            result.status = "partial" if result.steps_completed else "failed"
        else:
            result.status = "success"

        logger.info(f"Pipeline complete: {result.status}")
        return result


# ═══════════════════════════════════════════════════════════════
#  Protein Analysis Pipeline
# ═══════════════════════════════════════════════════════════════

class ProteinAnalysisPipeline:
    """
    Protein property prediction pipeline using XtriMoPGLM.

    Use Case:
    - Predict protein function, structure, and properties
    - Analyze amino acid sequences before simulation
    - Prepare data for EvoSim

    Examples:
        >>> pipeline = ProteinAnalysisPipeline(device="cuda")
        >>> result = pipeline.run("proteins.fasta", output_dir="/results")
        >>> print(f"Analyzed {len(result.results['predictions'])} proteins")
    """

    def __init__(
        self,
        model_path: Optional[str | Path] = None,
        device: str = "cuda",
        batch_size: int = 32,
        use_docker: bool = False,  # XtriMoPGLM typically runs locally with GPU
    ):
        """
        Initialize protein analysis pipeline.

        Args:
            model_path: Path to XtriMoPGLM model
            device: 'cuda' or 'cpu'
            batch_size: Batch size for inference
            use_docker: Use Docker (not recommended for GPU workloads)
        """
        config = XtriMoPGLMConfig(
            model_path=Path(model_path) if model_path else None,
            device=device,
            batch_size=batch_size,
            use_docker=use_docker,
        )
        self.xtrimopglm = XtriMoPGLMTool(config)

    def run(
        self,
        protein_fasta: str | Path,
        output_dir: Optional[str | Path] = None,
    ) -> PipelineResult:
        """
        Run protein analysis pipeline.

        Args:
            protein_fasta: FASTA file with protein sequences
            output_dir: Output directory

        Returns:
            PipelineResult with predictions
        """
        result = PipelineResult(
            pipeline_name="ProteinAnalysisPipeline",
            status="running"
        )

        logger.info("Starting protein analysis pipeline")

        try:
            logger.info("Predicting protein properties with XtriMoPGLM")
            prediction_result = self.xtrimopglm.predict(
                protein_fasta,
                output_dir=output_dir
            )

            result.steps_completed.append("xtrimopglm_prediction")
            result.results["predictions"] = prediction_result["results"]
            result.status = "success"

            logger.info("Protein analysis complete")

        except Exception as exc:
            logger.error(f"Protein analysis failed: {exc}")
            result.steps_failed.append("xtrimopglm_prediction")
            result.errors["xtrimopglm"] = str(exc)
            result.status = "failed"

        return result


# ═══════════════════════════════════════════════════════════════
#  Mutation Effect Prediction Pipeline
# ═══════════════════════════════════════════════════════════════

class MutationEffectPipeline:
    """
    Mutation effect prediction pipeline using Enformer/DeepSEA.

    Use Case:
    - Predict phenotypic effects of DNA mutations
    - Forecast gene expression changes
    - Integration with EvoSim for evolutionary simulations

    Examples:
        >>> pipeline = MutationEffectPipeline(model="enformer")
        >>> result = pipeline.run("mutations.fasta", output_dir="/results")
        >>> print(f"Predicted effects for {len(result.results['predictions'])} mutations")
    """

    def __init__(
        self,
        model_name: str = "enformer",  # 'enformer' or 'deepsea'
        device: str = "cuda",
        use_docker: bool = False,
    ):
        """
        Initialize mutation effect pipeline.

        Args:
            model_name: 'enformer' or 'deepsea'
            device: 'cuda' or 'cpu'
            use_docker: Use Docker
        """
        config = EnformerConfig(
            model_name=model_name,
            device=device,
            use_docker=use_docker,
        )
        self.enformer = EnformerTool(config)

    def run(
        self,
        sequence_fasta: str | Path,
        output_dir: Optional[str | Path] = None,
    ) -> PipelineResult:
        """
        Run mutation effect prediction pipeline.

        Args:
            sequence_fasta: FASTA file with DNA sequences (mutations)
            output_dir: Output directory

        Returns:
            PipelineResult with effect predictions
        """
        result = PipelineResult(
            pipeline_name="MutationEffectPipeline",
            status="running"
        )

        logger.info("Starting mutation effect prediction pipeline")

        try:
            logger.info("Predicting mutation effects with Enformer")
            prediction_result = self.enformer.predict(
                sequence_fasta,
                output_dir=output_dir
            )

            result.steps_completed.append("enformer_prediction")
            result.results["predictions"] = prediction_result["results"]
            result.status = "success"

            logger.info("Mutation effect prediction complete")

        except Exception as exc:
            logger.error(f"Mutation prediction failed: {exc}")
            result.steps_failed.append("enformer_prediction")
            result.errors["enformer"] = str(exc)
            result.status = "failed"

        return result


# ═══════════════════════════════════════════════════════════════
#  MetaGraph Data Source Integration
# ═══════════════════════════════════════════════════════════════

class MetaGraphDataSource:
    """
    MetaGraph integration as data source for GeneScope.

    Use Case:
    - Search for sequences across global databases
    - Find homologous sequences
    - Discover related genomic data

    Examples:
        >>> data_source = MetaGraphDataSource(
        ...     index_dir="/data/metagraph_sra_index"
        ... )
        >>> results = data_source.search_sequences(
        ...     "query.fasta",
        ...     discovery_fraction=0.8
        ... )
        >>> print(f"Found {results.results['num_matches']} matches")
    """

    def __init__(
        self,
        index_dir: str | Path,
        use_docker: bool = True,
    ):
        """
        Initialize MetaGraph data source.

        Args:
            index_dir: Path to MetaGraph index
            use_docker: Use Docker
        """
        config = MetaGraphConfig(
            index_dir=Path(index_dir),
            use_docker=use_docker,
        )
        self.metagraph = MetaGraphTool(config)

    def search_sequences(
        self,
        query_fasta: str | Path,
        output_dir: Optional[str | Path] = None,
        discovery_fraction: float = 0.8,
    ) -> PipelineResult:
        """
        Search for sequences in MetaGraph index.

        Args:
            query_fasta: FASTA file with query sequences
            output_dir: Output directory
            discovery_fraction: Fraction of k-mers to discover

        Returns:
            PipelineResult with search results
        """
        result = PipelineResult(
            pipeline_name="MetaGraphDataSource",
            status="running"
        )

        logger.info("Searching MetaGraph index")

        try:
            search_result = self.metagraph.search(
                query_fasta,
                output_dir=output_dir,
                discovery_fraction=discovery_fraction
            )

            result.steps_completed.append("metagraph_search")
            result.results.update(search_result["results"])
            result.status = "success"

            logger.info(f"MetaGraph search complete: {result.results['num_matches']} matches")

        except Exception as exc:
            logger.error(f"MetaGraph search failed: {exc}")
            result.steps_failed.append("metagraph_search")
            result.errors["metagraph"] = str(exc)
            result.status = "failed"

        return result

    def build_index(
        self,
        input_fasta: str | Path,
        output_dir: Optional[str | Path] = None,
        kmer_length: int = 31,
    ) -> PipelineResult:
        """
        Build MetaGraph index from sequences.

        Args:
            input_fasta: FASTA file with sequences
            output_dir: Output directory for index
            kmer_length: K-mer length

        Returns:
            PipelineResult with build status
        """
        result = PipelineResult(
            pipeline_name="MetaGraphIndexBuild",
            status="running"
        )

        logger.info("Building MetaGraph index")

        try:
            build_result = self.metagraph.build_index(
                input_fasta,
                output_dir=output_dir,
                kmer_length=kmer_length
            )

            result.steps_completed.append("metagraph_build")
            result.results.update(build_result["results"])
            result.status = "success"

            logger.info("MetaGraph index build complete")

        except Exception as exc:
            logger.error(f"MetaGraph build failed: {exc}")
            result.steps_failed.append("metagraph_build")
            result.errors["metagraph"] = str(exc)
            result.status = "failed"

        return result
