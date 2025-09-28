"""
Enhanced settings configuration for GenoScope with validation and environment support.
"""

from pathlib import Path

from pydantic import ConfigDict
from pydantic import Field
from pydantic import field_validator
from pydantic_settings import BaseSettings


class DatabaseSettings(BaseSettings):
    """Database configuration settings."""

    # SQLite settings
    database_url: str = Field(
        default="sqlite:///./genoscope.db",
        description="Database URL for SQLite or PostgreSQL",
    )

    # Connection settings
    pool_size: int = Field(default=5, ge=1, le=50)
    max_overflow: int = Field(default=10, ge=0, le=100)
    pool_timeout: int = Field(default=30, ge=1, le=300)
    pool_recycle: int = Field(default=3600, ge=300, le=86400)

    # SQLite specific
    sqlite_timeout: int = Field(default=20, ge=5, le=300)
    sqlite_check_same_thread: bool = False

    model_config = ConfigDict(env_prefix="DB_")


class SecuritySettings(BaseSettings):
    """Security configuration settings."""

    # File upload security
    max_file_size: int = Field(
        default=100 * 1024 * 1024, description="Maximum file size in bytes"  # 100MB
    )

    allowed_extensions: set[str] = Field(
        default={
            ".csv",
            ".tsv",
            ".vcf",
            ".parquet",
            ".bam",
            ".xlsx",
            ".xls",
            ".json",
            ".fasta",
            ".gff",
            ".hdf5",
        },
        description="Allowed file extensions",
    )

    # API security
    cors_origins: list[str] = Field(
        default=[
            "http://localhost:3000",
            "http://localhost:5173",
            "http://127.0.0.1:3000",
        ],
        description="Allowed CORS origins",
    )

    cors_allow_credentials: bool = False
    cors_allow_methods: list[str] = ["GET", "POST", "PUT", "DELETE", "OPTIONS"]
    cors_allow_headers: list[str] = ["*"]

    # Rate limiting
    rate_limit_enabled: bool = Field(
        default=True, description="Enable API rate limiting"
    )
    rate_limit_requests_per_minute: int = Field(default=60, ge=1, le=10000)
    rate_limit_burst: int = Field(default=100, ge=1, le=10000)

    # Authentication (for future use)
    secret_key: str | None = Field(default=None, description="JWT secret key")
    access_token_expire_minutes: int = Field(default=30, ge=5, le=1440)

    @field_validator("allowed_extensions", mode="before")
    @classmethod
    def normalize_extensions(cls, v):
        """Ensure all extensions start with a dot and are lowercase."""
        if isinstance(v, str):
            v = v.split(",")

        normalized = set()
        for ext in v:
            ext = ext.strip().lower()
            if not ext.startswith("."):
                ext = f".{ext}"
            normalized.add(ext)
        return normalized

    @field_validator("max_file_size")
    @classmethod
    def validate_file_size(cls, v):
        """Validate file size is reasonable."""
        min_size = 1024  # 1KB
        max_size = 1024 * 1024 * 1024 * 10  # 10GB

        if v < min_size:
            raise ValueError(f"max_file_size must be at least {min_size} bytes")
        if v > max_size:
            raise ValueError(f"max_file_size must not exceed {max_size} bytes")
        return v

    model_config = ConfigDict(env_prefix="SECURITY_")


class CacheSettings(BaseSettings):
    """Cache configuration settings."""

    # Redis settings
    redis_url: str = Field(default="redis://localhost:6379/0", description="Redis URL")
    redis_password: str | None = None
    redis_ssl: bool = False
    redis_timeout: int = Field(default=5, ge=1, le=60)

    # Cache behavior
    cache_enabled: bool = Field(default=True, description="Enable caching")
    default_cache_ttl: int = Field(default=3600, ge=60, le=86400)  # 1 hour default
    max_cache_size: int = Field(default=1000, ge=10, le=100000)

    # Specific cache TTLs
    analysis_result_ttl: int = Field(default=7200, ge=300, le=86400)  # 2 hours
    file_metadata_ttl: int = Field(default=1800, ge=60, le=3600)  # 30 minutes
    user_session_ttl: int = Field(default=1800, ge=300, le=86400)  # 30 minutes

    model_config = ConfigDict(env_prefix="CACHE_")


class CelerySettings(BaseSettings):
    """Celery task queue configuration."""

    # Broker settings
    broker_url: str = Field(
        default="redis://localhost:6379/1", description="Celery broker URL"
    )
    result_backend: str = Field(
        default="redis://localhost:6379/2", description="Celery result backend"
    )

    # Task settings
    task_serializer: str = "json"
    result_serializer: str = "json"
    accept_content: list[str] = ["json"]
    timezone: str = "UTC"
    enable_utc: bool = True

    # Worker settings
    worker_prefetch_multiplier: int = Field(default=1, ge=1, le=10)
    task_acks_late: bool = True
    task_reject_on_worker_lost: bool = True

    # Task routing
    task_routes: dict[str, dict[str, str]] = Field(
        default={
            "genoscope.tasks.analyze_file": {"queue": "analysis"},
            "genoscope.tasks.generate_report": {"queue": "reports"},
            "genoscope.tasks.cleanup_files": {"queue": "maintenance"},
        }
    )

    # Task time limits
    task_soft_time_limit: int = Field(default=3600, ge=60, le=86400)  # 1 hour
    task_time_limit: int = Field(default=4200, ge=120, le=86400)  # 1 hour 10 minutes

    model_config = ConfigDict(env_prefix="CELERY_")


class LoggingSettings(BaseSettings):
    """Logging configuration settings."""

    # Basic logging
    log_level: str = Field(default="INFO", description="Logging level")
    log_format: str = Field(
        default="detailed", description="Log format (simple, detailed, json)"
    )

    # File logging
    log_file: str | None = Field(default=None, description="Log file path")
    log_rotation_size: str = Field(
        default="100MB", description="Log file rotation size"
    )
    log_retention: str = Field(default="30 days", description="Log retention period")
    log_compression: bool = Field(default=True, description="Compress rotated logs")

    # Special loggers
    access_log_enabled: bool = Field(default=True, description="Enable access logging")
    performance_log_enabled: bool = Field(
        default=False, description="Enable performance logging"
    )
    sql_log_enabled: bool = Field(default=False, description="Enable SQL query logging")

    # External service logging
    uvicorn_log_level: str = Field(default="INFO", description="Uvicorn log level")
    sqlalchemy_log_level: str = Field(
        default="WARNING", description="SQLAlchemy log level"
    )

    @field_validator("log_level", "uvicorn_log_level", "sqlalchemy_log_level")
    @classmethod
    def validate_log_level(cls, v):
        """Validate log level is valid."""
        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        v = v.upper()
        if v not in valid_levels:
            raise ValueError(f"log_level must be one of {valid_levels}")
        return v

    @field_validator("log_format")
    @classmethod
    def validate_log_format(cls, v):
        """Validate log format is valid."""
        valid_formats = ["simple", "detailed", "json"]
        if v not in valid_formats:
            raise ValueError(f"log_format must be one of {valid_formats}")
        return v

    model_config = ConfigDict(env_prefix="LOG_")


class AnalysisSettings(BaseSettings):
    """Data analysis configuration settings."""

    # Processing limits
    max_rows_in_memory: int = Field(default=1000000, ge=1000, le=100000000)
    chunk_size: int = Field(default=10000, ge=100, le=1000000)
    max_file_processing_time: int = Field(default=3600, ge=60, le=86400)  # 1 hour

    # Analysis parameters
    default_pca_components: int = Field(default=2, ge=2, le=50)
    outlier_detection_method: str = Field(
        default="iqr", description="Default outlier detection method"
    )
    missing_value_method: str = Field(
        default="mean", description="Default missing value handling"
    )

    # Performance settings
    use_multiprocessing: bool = Field(
        default=True, description="Enable multiprocessing"
    )
    max_workers: int = Field(default=4, ge=1, le=32)
    memory_limit_gb: float = Field(default=4.0, ge=0.5, le=64.0)

    # Feature flags
    enable_ml_imputation: bool = Field(
        default=False, description="Enable ML-based imputation"
    )
    enable_advanced_stats: bool = Field(
        default=True, description="Enable advanced statistics"
    )
    enable_clustering: bool = Field(
        default=True, description="Enable clustering analysis"
    )

    @field_validator("outlier_detection_method")
    @classmethod
    def validate_outlier_method(cls, v):
        """Validate outlier detection method."""
        valid_methods = ["iqr", "z-score", "mahalanobis", "isolation_forest"]
        if v not in valid_methods:
            raise ValueError(f"outlier_detection_method must be one of {valid_methods}")
        return v

    @field_validator("missing_value_method")
    @classmethod
    def validate_missing_value_method(cls, v):
        """Validate missing value handling method."""
        valid_methods = [
            "mean",
            "median",
            "mode",
            "ffill",
            "bfill",
            "interpolate",
            "ml",
        ]
        if v not in valid_methods:
            raise ValueError(f"missing_value_method must be one of {valid_methods}")
        return v

    model_config = ConfigDict(env_prefix="ANALYSIS_")


class APISettings(BaseSettings):
    """API server configuration settings."""

    # Server settings
    host: str = Field(default="0.0.0.0", description="API server host")
    port: int = Field(default=8000, ge=1, le=65535)
    debug: bool = Field(default=False, description="Enable debug mode")
    reload: bool = Field(default=False, description="Enable auto-reload")

    # API behavior
    api_version: str = Field(default="v1", description="API version prefix")
    api_title: str = Field(default="GenoScope API", description="API title")
    api_description: str = Field(
        default="Advanced genomic data analysis platform", description="API description"
    )

    # Documentation
    docs_url: str = Field(default="/docs", description="Swagger UI docs URL")
    redoc_url: str = Field(default="/redoc", description="ReDoc docs URL")
    openapi_url: str = Field(default="/openapi.json", description="OpenAPI schema URL")

    # Request handling
    max_request_size: int = Field(
        default=100 * 1024 * 1024, ge=1024, le=1024 * 1024 * 1024
    )  # 100MB
    request_timeout: int = Field(default=300, ge=30, le=3600)  # 5 minutes

    model_config = ConfigDict(env_prefix="API_")


class DirectorySettings(BaseSettings):
    """Directory paths configuration."""

    # Base directories
    upload_dir: str = Field(default="uploads", description="Upload directory")
    data_dir: str = Field(default="data", description="Data directory")
    temp_dir: str = Field(default="temp", description="Temporary files directory")
    logs_dir: str = Field(default="logs", description="Logs directory")

    # Subdirectories
    datasets_dir: str = Field(
        default="data/datasets", description="Datasets storage directory"
    )
    reports_dir: str = Field(default="data/reports", description="Reports directory")
    cache_dir: str = Field(default="data/cache", description="Cache directory")
    backup_dir: str = Field(default="data/backups", description="Backup directory")

    # File cleanup
    cleanup_temp_files: bool = Field(
        default=True, description="Auto-cleanup temporary files"
    )
    temp_file_max_age_hours: int = Field(default=24, ge=1, le=168)  # 1 day to 1 week

    def ensure_directories(self) -> None:
        """Create all configured directories if they don't exist."""
        dirs_to_create = [
            self.upload_dir,
            self.data_dir,
            self.temp_dir,
            self.logs_dir,
            self.datasets_dir,
            self.reports_dir,
            self.cache_dir,
            self.backup_dir,
        ]

        for dir_path in dirs_to_create:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

    model_config = ConfigDict(env_prefix="DIR_")


class Settings(BaseSettings):
    """Main application settings."""

    # Environment
    env: str = Field(
        default="development",
        description="Environment (development, testing, production)",
    )
    debug: bool = Field(default=True, description="Debug mode")
    testing: bool = Field(default=False, description="Testing mode")

    # Subsettings
    database: DatabaseSettings = Field(default_factory=DatabaseSettings)
    security: SecuritySettings = Field(default_factory=SecuritySettings)
    cache: CacheSettings = Field(default_factory=CacheSettings)
    celery: CelerySettings = Field(default_factory=CelerySettings)
    logging: LoggingSettings = Field(default_factory=LoggingSettings)
    analysis: AnalysisSettings = Field(default_factory=AnalysisSettings)
    api: APISettings = Field(default_factory=APISettings)
    directories: DirectorySettings = Field(default_factory=DirectorySettings)

    # Legacy compatibility (for backward compatibility with existing code)
    @property
    def ENV(self) -> str:
        """Legacy compatibility property."""
        return self.env.upper()

    @property
    def UPLOAD_DIR(self) -> str:
        """Legacy compatibility property."""
        return self.directories.upload_dir

    @property
    def MAX_FILE_SIZE(self) -> int:
        """Legacy compatibility property."""
        return self.security.max_file_size

    @property
    def ALLOWED_EXTENSIONS(self) -> list[str]:
        """Legacy compatibility property."""
        return list(self.security.allowed_extensions)

    @property
    def ALLOW_ORIGINS(self) -> list[str]:
        """Legacy compatibility property."""
        return self.security.cors_origins

    @field_validator("env")
    @classmethod
    def validate_environment(cls, v):
        """Validate environment setting."""
        valid_envs = ["development", "testing", "staging", "production"]
        if v.lower() not in valid_envs:
            raise ValueError(f"env must be one of {valid_envs}")
        return v.lower()

    def is_development(self) -> bool:
        """Check if running in development mode."""
        return self.env == "development"

    def is_production(self) -> bool:
        """Check if running in production mode."""
        return self.env == "production"

    def is_testing(self) -> bool:
        """Check if running in testing mode."""
        return self.env == "testing" or self.testing

    def setup_environment(self) -> None:
        """Setup environment-specific configurations."""
        # Ensure directories exist
        self.directories.ensure_directories()

        # Adjust settings based on environment
        if self.is_production():
            self.debug = False
            self.api.debug = False
            self.api.reload = False
            self.logging.log_level = "INFO"
        elif self.is_development():
            self.debug = True
            self.api.debug = True
            self.api.reload = True
            self.logging.log_level = "DEBUG"
            self.logging.performance_log_enabled = True
        elif self.is_testing():
            self.debug = False
            self.logging.log_level = "WARNING"
            self.cache.cache_enabled = False

    model_config = ConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
        validate_assignment=True,
    )


# Global settings instance
settings = Settings()

# Setup environment on import
settings.setup_environment()


def get_settings() -> Settings:
    """Get the global settings instance."""
    return settings


def reload_settings() -> Settings:
    """Reload settings from environment."""
    global settings
    settings = Settings()
    settings.setup_environment()
    return settings


# Configuration validation
def validate_settings() -> None:
    """Validate all settings and raise errors for invalid configurations."""
    try:
        # Test database connection
        if settings.database.database_url.startswith("sqlite"):
            db_path = settings.database.database_url.replace("sqlite:///", "")
            if db_path != ":memory:":
                Path(db_path).parent.mkdir(parents=True, exist_ok=True)

        # Validate cache settings
        if settings.cache.cache_enabled and "redis" in settings.cache.redis_url:
            # Could add Redis connection test here
            pass

        # Validate directories
        settings.directories.ensure_directories()

        print("✅ Settings validation passed")

    except Exception as e:
        print(f"❌ Settings validation failed: {e}")
        raise


if __name__ == "__main__":
    # Validate settings when run directly
    validate_settings()

    # Print current configuration
    print("\n📋 Current Configuration:")
    print(f"Environment: {settings.env}")
    print(f"Debug mode: {settings.debug}")
    print(f"API host:port: {settings.api.host}:{settings.api.port}")
    print(f"Database: {settings.database.database_url}")
    print(f"Upload directory: {settings.directories.upload_dir}")
    print(f"Max file size: {settings.security.max_file_size / (1024*1024):.1f} MB")
    print(f"Cache enabled: {settings.cache.cache_enabled}")
    print(f"Log level: {settings.logging.log_level}")
