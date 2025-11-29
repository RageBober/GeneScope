"""
BioForge Configuration Settings.

Clean, validated configuration using pydantic-settings.
"""

from pathlib import Path
from typing import Annotated, Literal

from pydantic import BeforeValidator, ConfigDict, Field, field_validator
from pydantic_settings import BaseSettings


def parse_cors_origins(v):
    """Parse CORS origins from comma-separated string or list."""
    if isinstance(v, str):
        return [origin.strip() for origin in v.split(",") if origin.strip()]
    return v


CorsOriginsList = Annotated[list[str], BeforeValidator(parse_cors_origins)]


class DatabaseSettings(BaseSettings):
    """Database configuration."""
    
    url: str = Field(
        default="postgresql+asyncpg://bioforge:bioforge@localhost:5432/bioforge",
        description="Database URL (async)",
    )
    
    # Connection pool
    pool_size: int = Field(default=5, ge=1, le=50)
    max_overflow: int = Field(default=10, ge=0, le=100)
    pool_timeout: int = Field(default=30, ge=1, le=300)
    echo: bool = Field(default=False, description="Echo SQL queries")
    
    model_config = ConfigDict(env_prefix="DB_")


class RedisSettings(BaseSettings):
    """Redis configuration for cache and task queue."""
    
    url: str = Field(default="redis://localhost:6379/0")
    password: str | None = None
    ssl: bool = False
    timeout: int = Field(default=5, ge=1, le=60)
    
    # Cache settings
    cache_ttl: int = Field(default=3600, description="Default cache TTL in seconds")
    
    model_config = ConfigDict(env_prefix="REDIS_")


class StorageSettings(BaseSettings):
    """File storage configuration."""
    
    # Local directories
    upload_dir: Path = Field(default=Path("uploads"))
    data_dir: Path = Field(default=Path("data"))
    temp_dir: Path = Field(default=Path("temp"))
    
    # Subdirectories
    samples_dir: Path = Field(default=Path("data/samples"))
    results_dir: Path = Field(default=Path("data/results"))
    reports_dir: Path = Field(default=Path("data/reports"))
    references_dir: Path = Field(default=Path("data/references"))
    
    # File limits
    max_file_size: int = Field(
        default=500 * 1024 * 1024,  # 500MB
        description="Maximum upload file size in bytes"
    )
    
    allowed_extensions: set[str] = Field(
        default={
            ".fastq", ".fastq.gz", ".fq", ".fq.gz",
            ".bam", ".sam", ".cram",
            ".vcf", ".vcf.gz", ".bcf",
            ".fasta", ".fa", ".fna",
            ".gff", ".gff3", ".gtf",
            ".bed", ".bedgraph",
            ".csv", ".tsv", ".parquet",
        }
    )
    
    def ensure_directories(self) -> None:
        """Create all directories if they don't exist."""
        for attr in ['upload_dir', 'data_dir', 'temp_dir', 'samples_dir', 
                     'results_dir', 'reports_dir', 'references_dir']:
            path = getattr(self, attr)
            path.mkdir(parents=True, exist_ok=True)
    
    model_config = ConfigDict(env_prefix="STORAGE_")


class MinioSettings(BaseSettings):
    """MinIO/S3 object storage configuration."""
    
    endpoint: str = Field(default="localhost:9000")
    access_key: str = Field(default="bioforge")
    secret_key: str = Field(default="bioforge123")
    secure: bool = Field(default=False, description="Use HTTPS")
    
    # Buckets
    bucket_raw: str = Field(default="bioforge-raw")
    bucket_aligned: str = Field(default="bioforge-aligned")
    bucket_variants: str = Field(default="bioforge-variants")
    bucket_reports: str = Field(default="bioforge-reports")
    bucket_references: str = Field(default="bioforge-references")
    
    model_config = ConfigDict(env_prefix="MINIO_")


class SecuritySettings(BaseSettings):
    """Security configuration."""

    secret_key: str = Field(
        default="dev-secret-key-change-in-production",
        description="Secret key for JWT tokens"
    )
    algorithm: str = Field(default="HS256")
    access_token_expire_minutes: int = Field(default=30)
    refresh_token_expire_days: int = Field(default=7)

    # CORS
    cors_origins: CorsOriginsList = Field(
        default=["http://localhost:3000", "http://localhost:5173"]
    )
    cors_allow_credentials: bool = False

    # Rate limiting
    rate_limit_enabled: bool = True
    rate_limit_per_minute: int = Field(default=60)

    model_config = ConfigDict(env_prefix="SECURITY_")


class PipelineSettings(BaseSettings):
    """Bioinformatics pipeline configuration."""
    
    # Reference genome
    reference_genome: str = Field(default="GRCh38")
    reference_path: Path | None = None
    
    # Tool paths (Docker or local)
    use_docker: bool = Field(default=True)
    bwa_image: str = Field(default="biocontainers/bwa:v0.7.17_cv1")
    gatk_image: str = Field(default="broadinstitute/gatk:4.5.0.0")
    samtools_image: str = Field(default="biocontainers/samtools:v1.9-4-deb_cv1")
    
    # Processing
    threads: int = Field(default=4, ge=1, le=64)
    memory_gb: int = Field(default=8, ge=1, le=256)
    
    model_config = ConfigDict(env_prefix="PIPELINE_")


class Settings(BaseSettings):
    """Main application settings."""
    
    # Environment
    env: Literal["development", "testing", "staging", "production"] = Field(
        default="development"
    )
    debug: bool = Field(default=True)
    
    # App info
    app_name: str = Field(default="BioForge")
    app_version: str = Field(default="2.0.0")
    api_prefix: str = Field(default="/api/v1")
    
    # Server
    host: str = Field(default="0.0.0.0")
    port: int = Field(default=8000)
    
    # Subsettings
    database: DatabaseSettings = Field(default_factory=DatabaseSettings)
    redis: RedisSettings = Field(default_factory=RedisSettings)
    storage: StorageSettings = Field(default_factory=StorageSettings)
    minio: MinioSettings = Field(default_factory=MinioSettings)
    security: SecuritySettings = Field(default_factory=SecuritySettings)
    pipeline: PipelineSettings = Field(default_factory=PipelineSettings)
    
    # Legacy compatibility
    @property
    def ENV(self) -> str:
        return self.env.upper()
    
    @property
    def UPLOAD_DIR(self) -> str:
        return str(self.storage.upload_dir)
    
    @property
    def MAX_FILE_SIZE(self) -> int:
        return self.storage.max_file_size
    
    @property
    def ALLOWED_EXTENSIONS(self) -> list[str]:
        return list(self.storage.allowed_extensions)
    
    @property
    def ALLOW_ORIGINS(self) -> list[str]:
        return self.security.cors_origins
    
    def is_production(self) -> bool:
        return self.env == "production"
    
    def is_development(self) -> bool:
        return self.env == "development"
    
    def setup(self) -> None:
        """Setup environment."""
        self.storage.ensure_directories()
        
        if self.is_production():
            self.debug = False
    
    model_config = ConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )


# Global instance
try:
    settings = Settings()
    settings.setup()
except Exception as e:
    # Handle environment variable parsing errors
    import os
    # Clear problematic environment variables
    if 'SECURITY_CORS_ORIGINS' in os.environ:
        del os.environ['SECURITY_CORS_ORIGINS']
    if 'SECURITY_ALLOWED_EXTENSIONS' in os.environ:
        del os.environ['SECURITY_ALLOWED_EXTENSIONS']
    # Retry
    settings = Settings()
    settings.setup()


def get_settings() -> Settings:
    """Get settings instance (for dependency injection)."""
    return settings
