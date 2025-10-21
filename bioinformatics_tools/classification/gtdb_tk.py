"""
🦠 GTDB-Tk Integration

GTDB-Tk is a toolkit for classifying bacterial and archaeal genomes
using the Genome Database Taxonomy (GTDB).

Reference: https://github.com/Ecogenomics/GTDBTk
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

from ..base_tool import BaseBioTool
from ..config import GTDBTkConfig

logger = logging.getLogger(__name__)


class GTDBTkTool(BaseBioTool):
    """Wrapper for GTDB-Tk classifier"""

    def __init__(self, config: Optional[GTDBTkConfig] = None):
        if config is None:
            config = GTDBTkConfig()
        super().__init__(config)
        self.config: GTDBTkConfig = config

    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """Build GTDB-Tk classify command"""
        output_dir = output_path.parent / f"{input_path.stem}_gtdbtk"

        cmd = [
            "gtdbtk", "classify_wf",
            "--genome_dir", str(input_path.parent),
            "--out_dir", str(output_dir),
            "--cpus", str(self.config.num_threads),
            "--pplacer_cpus", str(self.config.pplacer_threads),
            "--min_af", str(self.config.min_af),
        ]

        if self.config.database_path:
            cmd.extend(["--db_path", str(self.config.database_path)])

        return cmd

    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """Parse GTDB-Tk output"""
        output_dir = output_path.parent / f"{output_path.stem}_gtdbtk"
        summary_file = output_dir / "gtdbtk.bac120.summary.tsv"

        if summary_file.exists():
            return {"classification_file": str(summary_file), "status": "success"}
        return {"status": "failed"}

    def classify(
        self,
        genome_dir: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Classify genomes using GTDB-Tk"""
        return self.execute(genome_dir, output_dir, **kwargs)
