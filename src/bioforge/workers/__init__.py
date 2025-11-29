"""ARQ Worker Configuration for BioForge.

Run worker with: uv run arq bioforge.workers.WorkerSettings
"""

from arq import cron
from arq.connections import RedisSettings

from bioforge.config import settings
from bioforge.workers.tasks import (
    cleanup_old_files,
    run_alignment,
    run_full_pipeline,
    run_qc,
    run_variant_calling,
)


def get_redis_settings() -> RedisSettings:
    """Get Redis settings from app config."""
    url = settings.redis.url

    # Parse redis://host:port/db URL
    if url.startswith("redis://"):
        url = url[8:]

    parts = url.split("/")
    host_port = parts[0]
    database = int(parts[1]) if len(parts) > 1 else 0

    if ":" in host_port:
        host, port_str = host_port.split(":")
        port = int(port_str)
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

    Available tasks:
    - run_qc: Quality control on FASTQ files
    - run_alignment: BWA-MEM2 alignment pipeline
    - run_variant_calling: bcftools/GATK variant calling
    - run_full_pipeline: Complete QC -> Alignment -> Variant Calling
    - cleanup_old_files: Remove old temp files
    """

    # Task functions available to workers
    functions = [
        run_qc,
        run_alignment,
        run_variant_calling,
        run_full_pipeline,
        cleanup_old_files,
    ]

    # Redis connection
    redis_settings = get_redis_settings()

    # Worker configuration
    max_jobs = 5  # Max concurrent jobs
    job_timeout = 7200  # 2 hours max per job
    keep_result = 86400  # Keep results for 24 hours
    keep_result_forever = False

    # Scheduled tasks (cron jobs)
    cron_jobs = [
        # Run cleanup at 3 AM daily
        cron(cleanup_old_files, hour=3, minute=0),
    ]

    # Retry configuration
    max_tries = 3
    retry_delay = 60  # 1 minute between retries

    # Health check
    health_check_interval = 60  # seconds

    @staticmethod
    async def on_startup(ctx: dict) -> None:
        """Called when worker starts."""
        from bioforge.core.logging import get_logger, setup_logging

        setup_logging(
            level="DEBUG" if settings.debug else "INFO",
            json_format=settings.is_production(),
        )

        logger = get_logger("bioforge.worker")
        logger.info(
            "BioForge ARQ Worker started",
            env=settings.env,
            max_jobs=WorkerSettings.max_jobs,
            job_timeout=WorkerSettings.job_timeout,
        )

        # Check available bioinformatics tools
        import shutil

        tools = {
            "bwa": shutil.which("bwa"),
            "bwa-mem2": shutil.which("bwa-mem2"),
            "samtools": shutil.which("samtools"),
            "bcftools": shutil.which("bcftools"),
            "gatk": shutil.which("gatk"),
            "fastp": shutil.which("fastp"),
        }

        for tool, path in tools.items():
            if path:
                logger.info(f"Tool available: {tool} -> {path}")
            else:
                logger.warning(f"Tool not found: {tool}")

        # Check if bioforge_core is available
        try:
            import bioforge_core

            logger.info(f"bioforge_core v{bioforge_core.__version__} available")
        except ImportError:
            logger.warning("bioforge_core not available (build with maturin)")

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


__all__ = [
    "WorkerSettings",
    "get_redis_settings",
    "run_qc",
    "run_alignment",
    "run_variant_calling",
    "run_full_pipeline",
    "cleanup_old_files",
]
