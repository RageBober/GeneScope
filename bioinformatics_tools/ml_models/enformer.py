"""
🧬 Enformer/DeepSEA Integration

Enformer is a deep learning model for predicting gene expression
and chromatin state from DNA sequence.

Reference: https://github.com/deepmind/deepmind-research/tree/master/enformer
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Optional

from ..base_tool import BaseBioTool
from ..config import EnformerConfig

logger = logging.getLogger(__name__)


class EnformerTool(BaseBioTool):
    """Wrapper for Enformer mutation effect prediction"""

    def __init__(self, config: Optional[EnformerConfig] = None):
        if config is None:
            config = EnformerConfig()
        super().__init__(config)
        self.config: EnformerConfig = config

    def _build_command(
        self,
        input_path: Path,
        output_path: Path,
        **kwargs: Any,
    ) -> list[str]:
        """Build Enformer prediction command"""
        cmd = [
            "python", "-m", "enformer.predict",
            "--input", str(input_path),
            "--output", str(output_path),
            "--model", self.config.model_name,
            "--device", self.config.device,
        ]

        return cmd

    def _parse_output(self, output_path: Path) -> dict[str, Any]:
        """Parse Enformer output"""
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
        """Predict mutation effects"""
        return self.execute(input_path, output_dir, **kwargs)
