"""
🤖 XtriMoPGLM Integration

XtriMoPGLM is a Transformer-based model for protein property prediction.

Use cases:
- Protein function prediction
- Structure prediction
- Sequence analysis

Reference: https://github.com/biomap-research/xTrimoGene
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Optional

from ..base_tool import BaseBioTool
from ..config import XtriMoPGLMConfig

logger = logging.getLogger(__name__)


class XtriMoPGLMTool(BaseBioTool):
    """Wrapper for XtriMoPGLM protein prediction model"""

    def __init__(self, config: Optional[XtriMoPGLMConfig] = None):
        if config is None:
            config = XtriMoPGLMConfig()
        super().__init__(config)
        self.config: XtriMoPGLMConfig = config

    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """Build XtriMoPGLM prediction command"""
        cmd = [
            "python", "-m", "xtrimopglm",
            "--input", str(input_path),
            "--output", str(output_path),
            "--batch_size", str(self.config.batch_size),
            "--device", self.config.device,
        ]

        if self.config.model_path:
            cmd.extend(["--model_path", str(self.config.model_path)])

        return cmd

    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """Parse XtriMoPGLM output"""
        if output_path.exists():
            try:
                with open(output_path) as f:
                    data = json.load(f)
                return {"predictions": data, "status": "success"}
            except Exception as exc:
                logger.error(f"Failed to parse output: {exc}")

        return {"status": "failed"}

    def predict(
        self,
        input_path: str | Path,
        output_dir: Optional[str | Path] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Predict protein properties"""
        return self.execute(input_path, output_dir, **kwargs)
