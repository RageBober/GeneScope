"""
⚙️ Configuration Management for Bioinformatics Tools

Centralized configuration for:
- Docker image versions
- Resource limits
- Database paths
- API endpoints
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class DockerConfig:
    """Docker container configuration"""

    # CPU limits
    cpu_count: int = 4
    cpu_quota: int = 100000  # microseconds per 100ms period

    # Memory limits
    memory_limit: str = "8g"
    memory_swap_limit: str = "16g"

    # Security
    read_only_root: bool = True
    no_new_privileges: bool = True
    user: str = "nobody"  # Run as non-root user

    # Network
    network_mode: str = "none"  # No network access by default

    # Timeout
    timeout_seconds: int = 3600  # 1 hour default


@dataclass
class ToolConfig:
    """Base configuration for bioinformatics tools"""

    # Tool metadata
    name: str
    version: str
    docker_image: Optional[str] = None

    # Resource limits
    docker_config: DockerConfig = field(default_factory=DockerConfig)

    # Paths
    work_dir: Path = field(default_factory=lambda: Path.cwd() / "work")
    cache_dir: Path = field(default_factory=lambda: Path.home() / ".cache/genoscope")

    # Execution mode
    use_docker: bool = True
    dry_run: bool = False

    def __post_init__(self):
        """Ensure directories exist"""
        self.work_dir = Path(self.work_dir)
        self.cache_dir = Path(self.cache_dir)

        for directory in [self.work_dir, self.cache_dir]:
            directory.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════
#  Tool-specific configurations
# ═══════════════════════════════════════════════════════════════

@dataclass
class MetaGraphConfig(ToolConfig):
    """Configuration for MetaGraph"""

    name: str = "MetaGraph"
    version: str = "0.3.0"
    docker_image: str = "quay.io/biocontainers/metagraph:0.3.0--py310h4de8cd1_0"

    # MetaGraph-specific settings
    index_dir: Optional[Path] = None
    kmer_length: int = 31
    bloom_filter_size: int = 10_000_000

    def __post_init__(self):
        super().__post_init__()
        if self.index_dir:
            self.index_dir = Path(self.index_dir)
            self.index_dir.mkdir(parents=True, exist_ok=True)


@dataclass
class Kraken2Config(ToolConfig):
    """Configuration for Kraken2/Centrifuge"""

    name: str = "Kraken2"
    version: str = "2.1.3"
    docker_image: str = "staphb/kraken2:2.1.3"

    # Kraken2-specific settings
    database_path: Optional[Path] = None
    confidence_threshold: float = 0.0
    quick_mode: bool = False
    min_hits: int = 1

    # Memory requirements for standard DB: ~50GB
    docker_config: DockerConfig = field(
        default_factory=lambda: DockerConfig(
            cpu_count=8,
            memory_limit="64g",
            memory_swap_limit="128g",
        )
    )

    def __post_init__(self):
        super().__post_init__()
        if self.database_path:
            self.database_path = Path(self.database_path)
            if not self.database_path.exists():
                raise ValueError(f"Kraken2 database not found: {self.database_path}")


@dataclass
class MEGAHITConfig(ToolConfig):
    """Configuration for MEGAHIT genome assembler"""

    name: str = "MEGAHIT"
    version: str = "1.2.9"
    docker_image: str = "vout/megahit:1.2.9"

    # MEGAHIT-specific settings
    min_contig_length: int = 200
    kmer_min: int = 21
    kmer_max: int = 141
    kmer_step: int = 12
    num_threads: int = 8

    # Memory requirements: depends on data size, ~100GB for large metagenomes
    docker_config: DockerConfig = field(
        default_factory=lambda: DockerConfig(
            cpu_count=16,
            memory_limit="128g",
            memory_swap_limit="256g",
            timeout_seconds=86400,  # 24 hours
        )
    )


@dataclass
class SPAdesConfig(ToolConfig):
    """Configuration for SPAdes genome assembler"""

    name: str = "SPAdes"
    version: str = "3.15.5"
    docker_image: str = "staphb/spades:3.15.5"

    # SPAdes-specific settings
    mode: str = "isolate"  # isolate, meta, rna, plasmid, etc.
    careful_mode: bool = True
    only_assembler: bool = False
    num_threads: int = 16

    docker_config: DockerConfig = field(
        default_factory=lambda: DockerConfig(
            cpu_count=16,
            memory_limit="128g",
            memory_swap_limit="256g",
            timeout_seconds=86400,  # 24 hours
        )
    )


@dataclass
class GTDBTkConfig(ToolConfig):
    """Configuration for GTDB-Tk bacterial classifier"""

    name: str = "GTDB-Tk"
    version: str = "2.3.2"
    docker_image: str = "ecogenomics/gtdbtk:2.3.2"

    # GTDB-Tk settings
    database_path: Optional[Path] = None
    num_threads: int = 8
    pplacer_threads: int = 4
    min_af: float = 0.65  # Minimum alignment fraction

    # GTDB database: ~80GB
    docker_config: DockerConfig = field(
        default_factory=lambda: DockerConfig(
            cpu_count=16,
            memory_limit="128g",
            memory_swap_limit="256g",
            timeout_seconds=14400,  # 4 hours
        )
    )

    def __post_init__(self):
        super().__post_init__()
        # GTDB database can be set via environment variable
        if not self.database_path:
            env_path = os.getenv("GTDBTK_DATA_PATH")
            if env_path:
                self.database_path = Path(env_path)

        if self.database_path:
            self.database_path = Path(self.database_path)
            if not self.database_path.exists():
                raise ValueError(f"GTDB-Tk database not found: {self.database_path}")


@dataclass
class XtriMoPGLMConfig(ToolConfig):
    """Configuration for XtriMoPGLM protein prediction model"""

    name: str = "XtriMoPGLM"
    version: str = "1.0.0"
    docker_image: Optional[str] = None  # Typically runs via Python/PyTorch

    # Model settings
    model_path: Optional[Path] = None
    device: str = "cuda"  # cuda or cpu
    batch_size: int = 32
    max_sequence_length: int = 1024

    # For GPU workloads
    docker_config: DockerConfig = field(
        default_factory=lambda: DockerConfig(
            cpu_count=8,
            memory_limit="32g",
            memory_swap_limit="64g",
            timeout_seconds=7200,  # 2 hours
        )
    )

    def __post_init__(self):
        super().__post_init__()
        if self.model_path:
            self.model_path = Path(self.model_path)
            if not self.model_path.exists():
                raise ValueError(f"XtriMoPGLM model not found: {self.model_path}")


@dataclass
class EnformerConfig(ToolConfig):
    """Configuration for Enformer/DeepSEA mutation effect prediction"""

    name: str = "Enformer"
    version: str = "1.0.0"
    docker_image: Optional[str] = None  # TensorFlow/JAX based

    # Model settings
    model_name: str = "enformer"  # enformer or deepsea
    sequence_length: int = 393216  # Enformer input length
    target_length: int = 896
    batch_size: int = 1
    device: str = "cuda"

    # Large memory requirements for genomic sequences
    docker_config: DockerConfig = field(
        default_factory=lambda: DockerConfig(
            cpu_count=16,
            memory_limit="64g",
            memory_swap_limit="128g",
            timeout_seconds=7200,  # 2 hours
        )
    )


# ═══════════════════════════════════════════════════════════════
#  Configuration Factory
# ═══════════════════════════════════════════════════════════════

def get_tool_config(tool_name: str, **kwargs) -> ToolConfig:
    """
    Factory function to get tool configuration by name.

    Args:
        tool_name: Name of the tool (case-insensitive)
        **kwargs: Override configuration parameters

    Returns:
        ToolConfig instance

    Raises:
        ValueError: If tool name is unknown

    Examples:
        >>> config = get_tool_config("kraken2", confidence_threshold=0.1)
        >>> config.name
        'Kraken2'
    """
    tool_configs = {
        "metagraph": MetaGraphConfig,
        "kraken2": Kraken2Config,
        "centrifuge": Kraken2Config,  # Similar config
        "megahit": MEGAHITConfig,
        "spades": SPAdesConfig,
        "gtdbtk": GTDBTkConfig,
        "gtdb-tk": GTDBTkConfig,
        "xtrimopglm": XtriMoPGLMConfig,
        "enformer": EnformerConfig,
        "deepsea": EnformerConfig,
    }

    tool_key = tool_name.lower().replace("_", "").replace("-", "")

    if tool_key not in tool_configs:
        raise ValueError(
            f"Unknown tool: {tool_name}. Available: {list(tool_configs.keys())}"
        )

    config_class = tool_configs[tool_key]
    return config_class(**kwargs)
