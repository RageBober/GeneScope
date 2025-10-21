"""
🦠 Kraken2 Integration

Kraken2 is a taxonomic classification system using exact k-mer matches
to achieve high accuracy and fast classification speeds.

Use cases:
- Taxonomic classification of metagenomic samples
- Contamination detection in genomic data
- Microbiome analysis

Reference: https://github.com/DerrickWood/kraken2
"""

from __future__ import annotations

import json
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Any, Optional

from ..base_tool import BaseBioTool, ToolExecutionError
from ..config import Kraken2Config

logger = logging.getLogger(__name__)


class Kraken2Tool(BaseBioTool):
    """
    Wrapper for Kraken2 taxonomic classification tool.

    Examples:
        >>> config = Kraken2Config(
        ...     database_path="/data/kraken2_db",
        ...     confidence_threshold=0.1
        ... )
        >>> tool = Kraken2Tool(config)
        >>> result = tool.classify("sample.fastq", output_dir="/results")
        >>> print(result['results']['classified_percentage'])
    """

    def __init__(self, config: Optional[Kraken2Config] = None):
        """
        Initialize Kraken2 tool.

        Args:
            config: Kraken2 configuration (default: Kraken2Config())
        """
        if config is None:
            config = Kraken2Config()
        super().__init__(config)
        self.config: Kraken2Config = config

    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """
        Build Kraken2 classification command.

        Args:
            input_path: Input sequences (FASTA/FASTQ format)
            output_path: Output file for classification results
            **kwargs: Additional parameters:
                - paired: Path to paired-end read file
                - report_file: Path for Kraken2 report
                - threads: Number of threads (default: from config)

        Returns:
            Command argument list
        """
        if not self.config.database_path:
            raise ValueError(
                "Kraken2 database path not set. "
                "Set via Kraken2Config(database_path='...') or KRAKEN2_DB_PATH env var"
            )

        # Kraken2 syntax:
        # kraken2 --db <db> --threads <n> --confidence <c> --output <out> --report <report> <input>

        report_file = kwargs.get("report_file", output_path.with_suffix(".report.txt"))
        threads = kwargs.get("threads", self.config.docker_config.cpu_count)

        cmd = [
            "kraken2",
            "--db", str(self.config.database_path),
            "--threads", str(threads),
            "--confidence", str(self.config.confidence_threshold),
            "--output", str(output_path),
            "--report", str(report_file),
        ]

        # Quick mode (less accurate but faster)
        if self.config.quick_mode:
            cmd.append("--quick")

        # Minimum hit groups
        if self.config.min_hits > 1:
            cmd.extend(["--minimum-hit-groups", str(self.config.min_hits)])

        # Paired-end reads
        paired_path = kwargs.get("paired")
        if paired_path:
            cmd.extend(["--paired", str(input_path), str(paired_path)])
        else:
            cmd.append(str(input_path))

        return cmd

    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """
        Parse Kraken2 output.

        Kraken2 outputs two files:
        1. Classification output (--output): per-read classification
        2. Report (--report): summary of taxa found

        Format of classification output:
        C/U <sequence_id> <taxon_id> <length> <kmer_info>

        Returns:
            Dictionary with parsed results
        """
        # Parse report file (more informative)
        report_path = output_path.with_suffix(".report.txt")

        results = {
            "classified_reads": 0,
            "unclassified_reads": 0,
            "total_reads": 0,
            "classified_percentage": 0.0,
            "taxa_found": [],
        }

        # Parse main output (classification per read)
        if output_path.exists():
            try:
                with open(output_path) as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue

                        parts = line.split("\t")
                        if len(parts) >= 1:
                            status = parts[0]
                            results["total_reads"] += 1

                            if status == "C":
                                results["classified_reads"] += 1
                            elif status == "U":
                                results["unclassified_reads"] += 1

            except Exception as exc:
                self.logger.error(f"Failed to parse Kraken2 output: {exc}")

        # Calculate percentage
        if results["total_reads"] > 0:
            results["classified_percentage"] = (
                100.0 * results["classified_reads"] / results["total_reads"]
            )

        # Parse report file (taxonomic breakdown)
        if report_path.exists():
            try:
                taxa = self._parse_report_file(report_path)
                results["taxa_found"] = taxa
            except Exception as exc:
                self.logger.error(f"Failed to parse Kraken2 report: {exc}")

        return results

    def _parse_report_file(self, report_path: Path) -> list[dict[str, Any]]:
        """
        Parse Kraken2 report file.

        Report format:
        <percentage> <num_clade_reads> <num_taxon_reads> <rank> <taxid> <name>

        Returns:
            List of taxa with their statistics
        """
        taxa = []

        with open(report_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Parse report line
                parts = line.split("\t")
                if len(parts) >= 6:
                    try:
                        percentage = float(parts[0])
                        num_clade_reads = int(parts[1])
                        num_taxon_reads = int(parts[2])
                        rank = parts[3].strip()
                        taxid = parts[4].strip()
                        name = parts[5].strip()

                        # Only include entries with actual reads
                        if num_taxon_reads > 0:
                            taxa.append({
                                "percentage": percentage,
                                "num_clade_reads": num_clade_reads,
                                "num_taxon_reads": num_taxon_reads,
                                "rank": rank,
                                "taxid": taxid,
                                "name": name,
                            })

                    except (ValueError, IndexError) as exc:
                        self.logger.warning(f"Skipping malformed report line: {line}")
                        continue

        # Sort by percentage (descending)
        taxa.sort(key=lambda x: x["percentage"], reverse=True)

        return taxa

    def classify(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Classify sequences using Kraken2.

        Args:
            input_path: Input sequences (FASTA/FASTQ)
            output_dir: Output directory
            **kwargs: Additional classification parameters

        Returns:
            Classification results dictionary
        """
        return self.execute(input_path, output_dir, **kwargs)


# ═══════════════════════════════════════════════════════════════
#  Convenience functions
# ═══════════════════════════════════════════════════════════════

def classify_kraken2(
    input_path: str | Path,
    database_path: str | Path,
    output_dir: Optional[str | Path] = None,
    confidence: float = 0.0,
    **kwargs: Any,
) -> dict[str, Any]:
    """
    Classify sequences with Kraken2 (convenience function).

    Args:
        input_path: Input sequences (FASTA/FASTQ)
        database_path: Kraken2 database directory
        output_dir: Output directory
        confidence: Confidence threshold (0.0-1.0)
        **kwargs: Additional parameters

    Returns:
        Classification results

    Examples:
        >>> results = classify_kraken2(
        ...     "sample.fastq",
        ...     database_path="/data/kraken2_std",
        ...     confidence=0.1
        ... )
        >>> print(f"Classified: {results['results']['classified_percentage']:.1f}%")
    """
    config = Kraken2Config(
        database_path=Path(database_path),
        confidence_threshold=confidence,
    )
    tool = Kraken2Tool(config)
    return tool.classify(input_path, output_dir, **kwargs)


def get_top_taxa(
    classification_results: dict[str, Any],
    n: int = 10,
    rank: Optional[str] = None,
) -> list[dict[str, Any]]:
    """
    Get top N classified taxa from Kraken2 results.

    Args:
        classification_results: Results from classify_kraken2()
        n: Number of top taxa to return
        rank: Filter by taxonomic rank (e.g., 'S' for species, 'G' for genus)

    Returns:
        List of top taxa

    Examples:
        >>> results = classify_kraken2("sample.fastq", "/data/kraken2_std")
        >>> top_species = get_top_taxa(results, n=5, rank='S')
        >>> for taxon in top_species:
        ...     print(f"{taxon['name']}: {taxon['percentage']:.2f}%")
    """
    if "results" not in classification_results:
        return []

    taxa = classification_results["results"].get("taxa_found", [])

    # Filter by rank if specified
    if rank:
        taxa = [t for t in taxa if t["rank"] == rank]

    # Return top N
    return taxa[:n]
