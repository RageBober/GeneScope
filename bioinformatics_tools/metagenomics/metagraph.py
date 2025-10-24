"""
🔍 MetaGraph Integration

MetaGraph is a tool for indexing and searching large collections of DNA/RNA sequences.
It builds compact searchable indexes of sequencing data.

Use cases:
- Search for sequences across multiple datasets
- Find similar sequences in large databases
- Metagenomic analysis

Reference: https://github.com/ratschlab/metagraph
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Optional

from ..base_tool import BaseBioTool, ToolExecutionError
from ..config import MetaGraphConfig

logger = logging.getLogger(__name__)


class MetaGraphTool(BaseBioTool):
    """
    Wrapper for MetaGraph sequence indexing and search tool.

    Examples:
        >>> config = MetaGraphConfig(index_dir="/data/metagraph_indices")
        >>> tool = MetaGraphTool(config)
        >>> result = tool.search("query.fasta", output_dir="/results")
        >>> print(result['results']['num_matches'])
    """

    def __init__(self, config: Optional[MetaGraphConfig] = None):
        """
        Initialize MetaGraph tool.

        Args:
            config: MetaGraph configuration (default: MetaGraphConfig())
        """
        if config is None:
            config = MetaGraphConfig()
        super().__init__(config)
        self.config: MetaGraphConfig = config

    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """
        Build MetaGraph search command.

        Args:
            input_path: Query sequences (FASTA format)
            output_path: Output file for matches
            **kwargs: Additional parameters:
                - mode: 'search' or 'build' (default: search)
                - discovery_fraction: Fraction of k-mers to discover (default: 0.8)
                - min_kmers_fraction: Min fraction of k-mers to match (default: 0.5)

        Returns:
            Command argument list
        """
        mode = kwargs.get("mode", "search")

        if mode == "search":
            return self._build_search_command(input_path, output_path, **kwargs)
        elif mode == "build":
            return self._build_index_command(input_path, output_path, **kwargs)
        else:
            raise ValueError(f"Unknown mode: {mode}. Use 'search' or 'build'")

    def _build_search_command(
        self,
        query_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """Build search command"""
        # MetaGraph search syntax:
        # metagraph query -i <index> --discovery-fraction <f> <query.fasta> > output.txt

        if not self.config.index_dir or not self.config.index_dir.exists():
            raise ValueError(
                f"MetaGraph index not found: {self.config.index_dir}. "
                "Run with mode='build' first."
            )

        discovery_fraction = kwargs.get("discovery_fraction", 0.8)
        min_kmers_fraction = kwargs.get("min_kmers_fraction", 0.5)

        cmd = [
            "metagraph", "query",
            "--index", str(self.config.index_dir / "graph"),
            "--discovery-fraction", str(discovery_fraction),
            "--min-kmers-fraction-label", str(min_kmers_fraction),
            "--output", str(output_path),
            str(query_path),
        ]

        return cmd

    def _build_index_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """Build indexing command"""
        # MetaGraph build syntax:
        # metagraph build -k <kmer_len> -o <output_prefix> <input.fasta>

        kmer_length = kwargs.get("kmer_length", self.config.kmer_length)

        if not self.config.index_dir:
            self.config.index_dir = output_path.parent / "metagraph_index"
            self.config.index_dir.mkdir(parents=True, exist_ok=True)

        output_prefix = self.config.index_dir / "graph"

        cmd = [
            "metagraph", "build",
            "-k", str(kmer_length),
            "--mode", "canonical",  # Use canonical k-mers
            "-o", str(output_prefix),
            str(input_path),
        ]

        return cmd

    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """
        Parse MetaGraph output.

        MetaGraph outputs tab-separated values with:
        - Query sequence ID
        - Matched labels/datasets
        - Match statistics

        Returns:
            Dictionary with parsed results
        """
        if not output_path.exists():
            self.logger.warning(f"Output file not found: {output_path}")
            return {"matches": [], "num_matches": 0}

        matches = []

        try:
            with open(output_path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    parts = line.split("\t")
                    if len(parts) >= 2:
                        matches.append({
                            "query_id": parts[0],
                            "matched_labels": parts[1].split(",") if parts[1] else [],
                            "raw_line": line,
                        })

        except Exception as exc:
            self.logger.error(f"Failed to parse MetaGraph output: {exc}")
            raise ToolExecutionError(f"Output parsing failed: {exc}") from exc

        return {
            "matches": matches,
            "num_matches": len(matches),
            "output_file": str(output_path),
        }

    def search(
        self,
        query_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Search for sequences in MetaGraph index.

        Args:
            query_path: Query sequences (FASTA format)
            output_dir: Output directory
            **kwargs: Additional search parameters

        Returns:
            Search results dictionary
        """
        kwargs["mode"] = "search"
        return self.execute(query_path, output_dir, **kwargs)

    def build_index(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Build MetaGraph index from input sequences.

        Args:
            input_path: Input sequences (FASTA format)
            output_dir: Output directory for index
            **kwargs: Additional build parameters

        Returns:
            Build results dictionary
        """
        kwargs["mode"] = "build"
        return self.execute(input_path, output_dir, **kwargs)


# ═══════════════════════════════════════════════════════════════
#  Convenience functions
# ═══════════════════════════════════════════════════════════════

def search_metagraph(
    query_path: str | Path,
    index_dir: str | Path,
    output_dir: Optional[str | Path] = None,
    **kwargs: Any,
) -> dict[str, Any]:
    """
    Search MetaGraph index (convenience function).

    Args:
        query_path: Query sequences (FASTA)
        index_dir: MetaGraph index directory
        output_dir: Output directory
        **kwargs: Additional parameters

    Returns:
        Search results

    Examples:
        >>> results = search_metagraph(
        ...     "queries.fasta",
        ...     index_dir="/data/metagraph_index",
        ...     discovery_fraction=0.9
        ... )
    """
    config = MetaGraphConfig(index_dir=Path(index_dir))
    tool = MetaGraphTool(config)
    return tool.search(query_path, output_dir, **kwargs)


def build_metagraph_index(
    input_path: str | Path,
    index_dir: str | Path,
    kmer_length: int = 31,
    **kwargs: Any,
) -> dict[str, Any]:
    """
    Build MetaGraph index (convenience function).

    Args:
        input_path: Input sequences (FASTA)
        index_dir: Output index directory
        kmer_length: K-mer length
        **kwargs: Additional parameters

    Returns:
        Build results

    Examples:
        >>> results = build_metagraph_index(
        ...     "genomes.fasta",
        ...     index_dir="/data/metagraph_index",
        ...     kmer_length=31
        ... )
    """
    config = MetaGraphConfig(index_dir=Path(index_dir), kmer_length=kmer_length)
    tool = MetaGraphTool(config)
    return tool.build_index(input_path, index_dir, **kwargs)
