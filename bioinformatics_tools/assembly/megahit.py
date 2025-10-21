"""
🧬 MEGAHIT Integration

MEGAHIT is an ultra-fast and memory-efficient NGS assembler optimized for
metagenomes, but also works well on generic single genome assembly.

Use cases:
- Metagenomic assembly from short reads
- Low-memory genome assembly
- Fast assembly of complex microbial communities

Reference: https://github.com/voutcn/megahit
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Optional

from ..base_tool import BaseBioTool, ToolExecutionError
from ..config import MEGAHITConfig

logger = logging.getLogger(__name__)


class MEGAHITTool(BaseBioTool):
    """
    Wrapper for MEGAHIT genome assembler.

    Examples:
        >>> config = MEGAHITConfig(min_contig_length=500, num_threads=16)
        >>> tool = MEGAHITTool(config)
        >>> result = tool.assemble("reads.fastq", output_dir="/results")
        >>> print(result['results']['num_contigs'])
    """

    def __init__(self, config: Optional[MEGAHITConfig] = None):
        """
        Initialize MEGAHIT tool.

        Args:
            config: MEGAHIT configuration (default: MEGAHITConfig())
        """
        if config is None:
            config = MEGAHITConfig()
        super().__init__(config)
        self.config: MEGAHITConfig = config

    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """
        Build MEGAHIT assembly command.

        Args:
            input_path: Input reads (FASTQ format)
            output_path: Output directory for assembly
            **kwargs: Additional parameters:
                - r2: Path to paired-end R2 reads
                - preset: Preset configuration (meta-sensitive, meta-large, etc.)

        Returns:
            Command argument list
        """
        # MEGAHIT syntax:
        # megahit -r <reads> -o <output_dir> --k-min <k> --k-max <k> --k-step <k> -t <threads>

        # Output directory (MEGAHIT requires a directory, not a file)
        output_dir = output_path.parent / f"{input_path.stem}_megahit_assembly"

        cmd = [
            "megahit",
            "--num-cpu-threads", str(self.config.num_threads),
            "--k-min", str(self.config.kmer_min),
            "--k-max", str(self.config.kmer_max),
            "--k-step", str(self.config.kmer_step),
            "--min-contig-len", str(self.config.min_contig_length),
            "--out-dir", str(output_dir),
        ]

        # Preset configuration
        preset = kwargs.get("preset")
        if preset:
            cmd.extend(["--presets", preset])

        # Input reads
        r2_path = kwargs.get("r2")
        if r2_path:
            # Paired-end
            cmd.extend([
                "-1", str(input_path),
                "-2", str(r2_path),
            ])
        else:
            # Single-end
            cmd.extend(["-r", str(input_path)])

        # Continue from previous run
        if kwargs.get("continue_run"):
            cmd.append("--continue")

        return cmd

    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """
        Parse MEGAHIT assembly output.

        MEGAHIT outputs:
        - final.contigs.fa: Assembled contigs
        - log: Assembly log with statistics

        Returns:
            Dictionary with parsed results
        """
        # Find output directory
        output_dir = output_path.parent / f"{output_path.stem}_megahit_assembly"

        results = {
            "contigs_file": None,
            "num_contigs": 0,
            "total_length": 0,
            "longest_contig": 0,
            "n50": 0,
        }

        # Check for contigs file
        contigs_file = output_dir / "final.contigs.fa"
        if contigs_file.exists():
            results["contigs_file"] = str(contigs_file)
            results.update(self._analyze_contigs(contigs_file))
        else:
            self.logger.warning(f"Contigs file not found: {contigs_file}")

        # Parse log file for additional statistics
        log_file = output_dir / "log"
        if log_file.exists():
            try:
                stats = self._parse_log(log_file)
                results.update(stats)
            except Exception as exc:
                self.logger.error(f"Failed to parse log: {exc}")

        return results

    def _analyze_contigs(self, contigs_file: Path) -> dict[str, Any]:
        """
        Analyze assembled contigs FASTA file.

        Calculates:
        - Number of contigs
        - Total assembly length
        - Longest contig
        - N50 statistic

        Returns:
            Dictionary with contig statistics
        """
        contig_lengths = []

        try:
            with open(contigs_file) as f:
                current_seq = ""

                for line in f:
                    line = line.strip()

                    if line.startswith(">"):
                        # New contig header
                        if current_seq:
                            contig_lengths.append(len(current_seq))
                            current_seq = ""
                    else:
                        # Sequence data
                        current_seq += line

                # Don't forget last contig
                if current_seq:
                    contig_lengths.append(len(current_seq))

        except Exception as exc:
            self.logger.error(f"Failed to analyze contigs: {exc}")
            return {}

        # Calculate statistics
        if not contig_lengths:
            return {
                "num_contigs": 0,
                "total_length": 0,
                "longest_contig": 0,
                "n50": 0,
            }

        contig_lengths.sort(reverse=True)

        total_length = sum(contig_lengths)
        half_length = total_length / 2

        # Calculate N50
        cumulative = 0
        n50 = 0
        for length in contig_lengths:
            cumulative += length
            if cumulative >= half_length:
                n50 = length
                break

        return {
            "num_contigs": len(contig_lengths),
            "total_length": total_length,
            "longest_contig": contig_lengths[0],
            "n50": n50,
        }

    def _parse_log(self, log_file: Path) -> dict[str, Any]:
        """
        Parse MEGAHIT log file for additional statistics.

        Returns:
            Dictionary with log statistics
        """
        stats = {}

        try:
            with open(log_file) as f:
                log_content = f.read()

                # Extract key statistics using regex
                # Example: "Total time: 123.45 seconds"
                time_match = re.search(r"Total time:\s+([\d.]+)", log_content)
                if time_match:
                    stats["assembly_time_seconds"] = float(time_match.group(1))

                # Example: "Peak memory usage: 12.34 GB"
                memory_match = re.search(r"Peak memory usage:\s+([\d.]+)\s+GB", log_content)
                if memory_match:
                    stats["peak_memory_gb"] = float(memory_match.group(1))

        except Exception as exc:
            self.logger.error(f"Failed to parse log: {exc}")

        return stats

    def assemble(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Assemble genome/metagenome using MEGAHIT.

        Args:
            input_path: Input reads (FASTQ)
            output_dir: Output directory
            **kwargs: Additional assembly parameters

        Returns:
            Assembly results dictionary
        """
        return self.execute(input_path, output_dir, **kwargs)


# ═══════════════════════════════════════════════════════════════
#  Convenience functions
# ═══════════════════════════════════════════════════════════════

def assemble_megahit(
    input_path: str | Path,
    output_dir: Optional[str | Path] = None,
    r2_path: Optional[str | Path] = None,
    preset: str = "meta-sensitive",
    **kwargs: Any,
) -> dict[str, Any]:
    """
    Assemble genome with MEGAHIT (convenience function).

    Args:
        input_path: Input reads (FASTQ)
        output_dir: Output directory
        r2_path: Paired-end R2 reads (optional)
        preset: MEGAHIT preset (meta-sensitive, meta-large, etc.)
        **kwargs: Additional parameters

    Returns:
        Assembly results

    Examples:
        >>> results = assemble_megahit(
        ...     "reads_R1.fastq",
        ...     r2_path="reads_R2.fastq",
        ...     preset="meta-sensitive"
        ... )
        >>> print(f"Assembled {results['results']['num_contigs']} contigs")
    """
    config = MEGAHITConfig()
    tool = MEGAHITTool(config)

    kwargs["r2"] = r2_path
    kwargs["preset"] = preset

    return tool.assemble(input_path, output_dir, **kwargs)
