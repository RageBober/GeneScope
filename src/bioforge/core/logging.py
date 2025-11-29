"""
Structured Logging Configuration.

Uses structlog for structured, JSON-formatted logging with OpenTelemetry integration.
"""

import logging
import sys
from typing import Any

import structlog
from structlog.types import Processor


def setup_logging(
    level: str = "INFO",
    json_format: bool = False,
    add_timestamp: bool = True,
) -> None:
    """
    Configure structlog for the application.
    
    Args:
        level: Log level (DEBUG, INFO, WARNING, ERROR)
        json_format: Use JSON output (for production)
        add_timestamp: Add timestamp to logs
    """
    # Shared processors
    shared_processors: list[Processor] = [
        structlog.contextvars.merge_contextvars,
        structlog.stdlib.add_log_level,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.StackInfoRenderer(),
        structlog.processors.UnicodeDecoder(),
    ]
    
    if add_timestamp:
        shared_processors.insert(0, structlog.processors.TimeStamper(fmt="iso"))
    
    if json_format:
        # Production: JSON output
        shared_processors.append(structlog.processors.format_exc_info)
        renderer: Processor = structlog.processors.JSONRenderer()
    else:
        # Development: colored console output
        renderer = structlog.dev.ConsoleRenderer(colors=True)
    
    structlog.configure(
        processors=shared_processors + [
            structlog.stdlib.ProcessorFormatter.wrap_for_formatter,
        ],
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )
    
    # Configure standard logging
    formatter = structlog.stdlib.ProcessorFormatter(
        foreign_pre_chain=shared_processors,
        processors=[
            structlog.stdlib.ProcessorFormatter.remove_processors_meta,
            renderer,
        ],
    )
    
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)
    
    # Root logger
    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    root_logger.addHandler(handler)
    root_logger.setLevel(getattr(logging, level.upper()))
    
    # Set levels for noisy loggers
    logging.getLogger("uvicorn.access").setLevel(logging.WARNING)
    logging.getLogger("sqlalchemy.engine").setLevel(logging.WARNING)
    logging.getLogger("httpx").setLevel(logging.WARNING)
    logging.getLogger("httpcore").setLevel(logging.WARNING)


def get_logger(name: str | None = None) -> structlog.stdlib.BoundLogger:
    """
    Get a structured logger.
    
    Args:
        name: Logger name (usually __name__)
    
    Returns:
        Configured structlog logger
    """
    return structlog.get_logger(name)


def bind_context(**kwargs: Any) -> None:
    """
    Bind context variables to all subsequent log messages.
    
    Useful for request-scoped data like request_id, user_id, etc.
    """
    structlog.contextvars.bind_contextvars(**kwargs)


def clear_context() -> None:
    """Clear all bound context variables."""
    structlog.contextvars.clear_contextvars()


# Request logging middleware helper
class RequestContextMiddleware:
    """Middleware to add request context to logs."""
    
    async def __call__(self, request, call_next):
        import uuid
        from starlette.requests import Request
        
        request_id = str(uuid.uuid4())[:8]
        
        bind_context(
            request_id=request_id,
            method=request.method,
            path=request.url.path,
        )
        
        try:
            response = await call_next(request)
            return response
        finally:
            clear_context()
