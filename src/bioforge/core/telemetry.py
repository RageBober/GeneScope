"""
OpenTelemetry Configuration.

Provides distributed tracing and metrics collection.
"""

from typing import Optional

from opentelemetry import trace
from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter
from opentelemetry.instrumentation.fastapi import FastAPIInstrumentor
from opentelemetry.instrumentation.sqlalchemy import SQLAlchemyInstrumentor
from opentelemetry.instrumentation.redis import RedisInstrumentor
from opentelemetry.sdk.resources import Resource, SERVICE_NAME, SERVICE_VERSION
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.trace.export import BatchSpanProcessor, ConsoleSpanExporter

from bioforge.config import settings
from bioforge.core.logging import get_logger

logger = get_logger(__name__)


def setup_telemetry(
    service_name: str = "bioforge",
    service_version: str = "2.0.0",
    otlp_endpoint: Optional[str] = None,
    console_export: bool = False,
) -> None:
    """
    Configure OpenTelemetry tracing.
    
    Args:
        service_name: Name of the service
        service_version: Version of the service
        otlp_endpoint: OTLP collector endpoint (e.g., "http://localhost:4317")
        console_export: Export traces to console (for development)
    """
    # Create resource
    resource = Resource.create({
        SERVICE_NAME: service_name,
        SERVICE_VERSION: service_version,
        "deployment.environment": settings.env,
    })
    
    # Create tracer provider
    provider = TracerProvider(resource=resource)
    
    # Add exporters
    if console_export or settings.debug:
        provider.add_span_processor(
            BatchSpanProcessor(ConsoleSpanExporter())
        )
    
    if otlp_endpoint:
        otlp_exporter = OTLPSpanExporter(endpoint=otlp_endpoint, insecure=True)
        provider.add_span_processor(BatchSpanProcessor(otlp_exporter))
        logger.info("OTLP exporter configured", endpoint=otlp_endpoint)
    
    # Set global tracer provider
    trace.set_tracer_provider(provider)
    
    logger.info(
        "Telemetry configured",
        service=service_name,
        version=service_version,
    )


def instrument_fastapi(app) -> None:
    """Instrument FastAPI application."""
    FastAPIInstrumentor.instrument_app(
        app,
        excluded_urls="health,healthz,ready,metrics",
    )
    logger.info("FastAPI instrumented")


def instrument_sqlalchemy(engine) -> None:
    """Instrument SQLAlchemy engine."""
    SQLAlchemyInstrumentor().instrument(
        engine=engine.sync_engine,
        enable_commenter=True,
    )
    logger.info("SQLAlchemy instrumented")


def instrument_redis() -> None:
    """Instrument Redis client."""
    RedisInstrumentor().instrument()
    logger.info("Redis instrumented")


def get_tracer(name: str = __name__) -> trace.Tracer:
    """Get a tracer instance."""
    return trace.get_tracer(name)


# Convenience decorator for tracing functions
def traced(name: Optional[str] = None):
    """
    Decorator to trace a function.
    
    Usage:
        @traced("my-operation")
        async def my_function():
            ...
    """
    def decorator(func):
        import functools
        
        span_name = name or func.__name__
        tracer = get_tracer(func.__module__)
        
        if asyncio_iscoroutinefunction(func):
            @functools.wraps(func)
            async def async_wrapper(*args, **kwargs):
                with tracer.start_as_current_span(span_name) as span:
                    span.set_attribute("function.name", func.__name__)
                    span.set_attribute("function.module", func.__module__)
                    try:
                        return await func(*args, **kwargs)
                    except Exception as e:
                        span.record_exception(e)
                        span.set_status(trace.Status(trace.StatusCode.ERROR))
                        raise
            return async_wrapper
        else:
            @functools.wraps(func)
            def sync_wrapper(*args, **kwargs):
                with tracer.start_as_current_span(span_name) as span:
                    span.set_attribute("function.name", func.__name__)
                    span.set_attribute("function.module", func.__module__)
                    try:
                        return func(*args, **kwargs)
                    except Exception as e:
                        span.record_exception(e)
                        span.set_status(trace.Status(trace.StatusCode.ERROR))
                        raise
            return sync_wrapper
    
    return decorator


def asyncio_iscoroutinefunction(func) -> bool:
    """Check if function is a coroutine."""
    import asyncio
    import inspect
    return asyncio.iscoroutinefunction(func) or inspect.iscoroutinefunction(func)
