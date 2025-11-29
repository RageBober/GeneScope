"""ARQ Worker runner script."""

import asyncio
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("bioforge.worker")


async def main():
    from arq import Worker
    from bioforge.workers import get_redis_settings
    from bioforge.workers.tasks import (
        run_alignment,
        run_variant_calling,
        run_full_pipeline,
        cleanup_old_files,
    )
    
    logger.info("Initializing worker...")
    
    redis_settings = get_redis_settings()
    logger.info(f"Redis: {redis_settings.host}:{redis_settings.port}")
    
    worker = Worker(
        functions=[
            run_alignment,
            run_variant_calling,
            run_full_pipeline,
            cleanup_old_files,
        ],
        redis_settings=redis_settings,
        max_jobs=10,
        job_timeout=3600,
        keep_result=3600,
    )
    
    logger.info("Worker starting...")
    await worker.main()


if __name__ == "__main__":
    print("=" * 50)
    print("BioForge ARQ Worker")
    print("=" * 50)
    asyncio.run(main())
