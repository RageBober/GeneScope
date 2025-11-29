"""ARQ Worker Configuration."""

from arq import cron
from arq.connections import RedisSettings

from bioforge.config import settings
from bioforge.workers.tasks import (
    run_alignment,
    run_variant_calling,
    run_full_pipeline,
    cleanup_old_files,
)


def get_redis_settings() -> RedisSettings:
    """Get Redis settings from app config."""
    # Parse redis URL
    url = settings.redis.url
    # redis://localhost:6379/0
    
    if url.startswith("redis://"):
        url = url[8:]
    
    parts = url.split("/")
    host_port = parts[0]
    database = int(parts[1]) if len(parts) > 1 else 0
    
    if ":" in host_port:
        host, port = host_port.split(":")
        port = int(port)
    else:
        host = host_port
        port = 6379
    
    return RedisSettings(
        host=host,
        port=port,
        database=database,
        password=settings.redis.password,
    )


class WorkerSettings:
    """ARQ Worker Settings.
    
    Run with: uv run arq bioforge.workers.WorkerSettings
    """
    
    # Task functions
    functions = [
        run_alignment,
        run_variant_calling,
        run_full_pipeline,
        cleanup_old_files,
    ]
    
    # Redis connection
    redis_settings = get_redis_settings()
    
    # Worker configuration
    max_jobs = 10
    job_timeout = 3600  # 1 hour
    keep_result = 3600  # Keep results for 1 hour
    
    # Scheduled tasks (cron jobs)
    cron_jobs = [
        cron(cleanup_old_files, hour=3, minute=0),  # Run cleanup at 3 AM
    ]
    
    # Retry configuration
    max_tries = 3
    retry_delay = 60  # 1 minute between retries
    
    # Health check
    health_check_interval = 60  # seconds
    
    @staticmethod
    async def on_startup(ctx: dict) -> None:
        """Called when worker starts."""
        from bioforge.core.logging import setup_logging, get_logger
        
        setup_logging(
            level="DEBUG" if settings.debug else "INFO",
            json_format=settings.is_production(),
        )
        
        logger = get_logger("bioforge.worker")
        logger.info(
            "BioForge ARQ Worker started",
            env=settings.env,
            max_jobs=WorkerSettings.max_jobs,
        )
    
    @staticmethod
    async def on_shutdown(ctx: dict) -> None:
        """Called when worker shuts down."""
        from bioforge.core.logging import get_logger
        logger = get_logger("bioforge.worker")
        logger.info("BioForge ARQ Worker shutting down")
    
    @staticmethod
    async def on_job_start(ctx: dict) -> None:
        """Called when a job starts."""
        from bioforge.core.logging import bind_context
        job_id = ctx.get("job_id", "unknown")
        bind_context(job_id=job_id)
    
    @staticmethod
    async def on_job_end(ctx: dict) -> None:
        """Called when a job ends."""
        from bioforge.core.logging import clear_context
        clear_context()
