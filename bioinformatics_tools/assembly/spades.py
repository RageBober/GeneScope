"""
🧬 SPAdes Integration

SPAdes is a versatile genome assembler for isolates, single-cell, metagenomes,
plasmids, RNA-Seq, and more.

Reference: https://github.com/ablab/spades
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

from ..base_tool import BaseBioTool
from ..config import SPAdesConfig

logger = logging.getLogger(__name__)


class SPAdesTool(BaseBioTool):
    """Wrapper for SPAdes genome assembler"""

    def __init__(self, config: Optional[SPAdesConfig] = None):
        if config is None:
            config = SPAdesConfig()
        super().__init__(config)
        self.config: SPAdesConfig = config

    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """Build SPAdes command"""
        output_dir = output_path.parent / f"{input_path.stem}_spades_assembly"

        cmd = [
            "spades.py",
            "-o", str(output_dir),
            "-t", str(self.config.num_threads),
        ]

        # Mode selection
        if self.config.mode == "meta":
            cmd.append("--meta")
        elif self.config.mode == "rna":
            cmd.append("--rna")
        elif self.config.mode == "plasmid":
            cmd.append("--plasmid")

        # Careful mode
        if self.config.careful_mode:
            cmd.append("--careful")

        # Input reads
        r2_path = kwargs.get("r2")
        if r2_path:
            cmd.extend(["-1", str(input_path), "-2", str(r2_path)])
        else:
            cmd.extend(["-s", str(input_path)])

        return cmd

    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """Parse SPAdes output"""
        output_dir = output_path.parent / f"{output_path.stem}_spades_assembly"
        contigs_file = output_dir / "contigs.fasta"

        if contigs_file.exists():
            return {"contigs_file": str(contigs_file), "status": "success"}
        return {"status": "failed", "contigs_file": None}

    def assemble(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Assemble genome using SPAdes"""
        return self.execute(input_path, output_dir, **kwargs)
