"""FastAPI application factory."""

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import ORJSONResponse

from bioforge.config import settings
from bioforge.core.exceptions import install_exception_handlers
from bioforge.core.logging import setup_logging, get_logger, RequestContextMiddleware
from bioforge.db import init_db, close_db
from bioforge.api.deps import close_arq_pool


# Initialize logging early
setup_logging(
    level="DEBUG" if settings.debug else "INFO",
    json_format=settings.is_production(),
)

logger = get_logger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan: startup and shutdown."""
    logger.info(
        "Starting BioForge",
        version=settings.app_version,
        env=settings.env,
    )
    
    # Startup
    await init_db()
    
    # Setup telemetry in production (before other initializations)
    if settings.is_production():
        from bioforge.core.telemetry import (
            setup_telemetry,
            instrument_fastapi,
            instrument_sqlalchemy,
            instrument_redis,
        )
        from bioforge.db import engine
        
        otlp_endpoint = os.getenv("OTLP_ENDPOINT")
        if otlp_endpoint:
            setup_telemetry(
                service_name=settings.app_name.lower(),
                service_version=settings.app_version,
                otlp_endpoint=otlp_endpoint,
            )
            instrument_fastapi(app)
            instrument_sqlalchemy(engine)
            instrument_redis()
    
    # Initialize S3 buckets if MinIO is configured
    if not settings.is_production() or os.getenv("MINIO_ENDPOINT"):
        try:
            from bioforge.services.s3_async import init_s3
            await init_s3()
            logger.info("S3 storage initialized")
        except Exception as e:
            logger.warning("S3 initialization skipped", error=str(e))
    
    logger.info("BioForge started successfully")
    
    yield
    
    # Shutdown
    logger.info("Shutting down BioForge")
    await close_arq_pool()
    await close_db()
    logger.info("BioForge shutdown complete")


def create_app() -> FastAPI:
    """Create and configure FastAPI application."""
    
    app = FastAPI(
        title=settings.app_name,
        version=settings.app_version,
        description="Genomic Analysis Platform",
        default_response_class=ORJSONResponse,
        docs_url="/docs" if settings.debug else None,
        redoc_url="/redoc" if settings.debug else None,
        openapi_url="/openapi.json" if settings.debug else None,
        lifespan=lifespan,
    )
    
    # Exception handlers
    install_exception_handlers(app, include_trace=settings.debug)
    
    # Middleware (order matters - last added is first executed)
    app.add_middleware(GZipMiddleware, minimum_size=1000)
    
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.security.cors_origins,
        allow_credentials=settings.security.cors_allow_credentials,
        allow_methods=["*"],
        allow_headers=["*"],
    )
    
    # Request context middleware for structured logging
    @app.middleware("http")
    async def add_request_context(request, call_next):
        import uuid
        from bioforge.core.logging import bind_context, clear_context
        
        request_id = request.headers.get("X-Request-ID", str(uuid.uuid4())[:8])
        
        bind_context(
            request_id=request_id,
            method=request.method,
            path=request.url.path,
        )
        
        try:
            response = await call_next(request)
            response.headers["X-Request-ID"] = request_id
            return response
        finally:
            clear_context()
    
    # Include routers
    from bioforge.api.routes import (
        health_router,
        auth_router,
        projects_router,
        samples_router,
        jobs_router,
        variants_router,
    )
    
    app.include_router(health_router)
    app.include_router(auth_router, prefix=settings.api_prefix)
    app.include_router(projects_router, prefix=settings.api_prefix)
    app.include_router(samples_router, prefix=settings.api_prefix)
    app.include_router(jobs_router, prefix=settings.api_prefix)
    app.include_router(variants_router, prefix=settings.api_prefix)
    
    return app


# Application instance
app = create_app()
