"""
Enhanced logging configuration for GenoScope.
Provides structured logging with different levels and formatters.
"""

import logging
import logging.config
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


def setup_logging(
    log_level: str = "INFO",
    log_file: Optional[str] = None,
    enable_structured: bool = False,
    enable_performance: bool = False,
) -> None:
    """
    Setup comprehensive logging configuration.

    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional log file path
        enable_structured: Enable structured JSON logging
        enable_performance: Enable performance timing logs
    """

    # Create logs directory if it doesn't exist
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

    # Base configuration
    config: Dict[str, Any] = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "standard": {
                "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
            "detailed": {
                "format": "%(asctime)s - %(name)s - %(levelname)s - %(module)s - %(funcName)s:%(lineno)d - %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
            "json": {
                "format": "%(asctime)s %(name)s %(levelname)s %(module)s %(funcName)s %(lineno)d %(message)s"
            },
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "level": log_level,
                "formatter": "standard",
                "stream": sys.stdout,
            },
            "console_error": {
                "class": "logging.StreamHandler",
                "level": "ERROR",
                "formatter": "detailed",
                "stream": sys.stderr,
            },
        },
        "loggers": {
            "genoscope": {
                "level": log_level,
                "handlers": ["console"],
                "propagate": False,
            },
            "genoscope.error": {
                "level": "ERROR",
                "handlers": ["console_error"],
                "propagate": False,
            },
            "uvicorn": {"level": "INFO", "handlers": ["console"], "propagate": False},
            "fastapi": {"level": "INFO", "handlers": ["console"], "propagate": False},
        },
        "root": {"level": log_level, "handlers": ["console"]},
    }

    # Add file handler if log_file is specified
    if log_file:
        config["handlers"]["file"] = {
            "class": "logging.handlers.RotatingFileHandler",
            "level": log_level,
            "formatter": "json" if enable_structured else "detailed",
            "filename": log_file,
            "maxBytes": 10485760,  # 10MB
            "backupCount": 5,
            "encoding": "utf8",
        }

        # Add file handler to all loggers
        for logger_config in config["loggers"].values():
            logger_config["handlers"].append("file")
        config["root"]["handlers"].append("file")

    # Add performance logger if enabled
    if enable_performance:
        config["loggers"]["genoscope.performance"] = {
            "level": "DEBUG",
            "handlers": ["console"],
            "propagate": False,
        }

    logging.config.dictConfig(config)


def get_logger(name: str = __name__) -> logging.Logger:
    """Get a configured logger instance."""
    return logging.getLogger(name)


def get_request_logger():
    """Get a logger for HTTP requests."""
    return logging.getLogger("genoscope.requests")


def get_analysis_logger():
    """Get a logger for data analysis operations."""
    return logging.getLogger("genoscope.analysis")


def get_performance_logger():
    """Get a logger for performance monitoring."""
    return logging.getLogger("genoscope.performance")


class PerformanceTimer:
    """Context manager for timing operations."""

    def __init__(self, operation_name: str, logger: logging.Logger | None = None):
        self.operation_name = operation_name
        self.logger = logger or get_performance_logger()
        self.start_time: Optional[float] = None

    def __enter__(self):
        self.start_time = datetime.now().timestamp()
        self.logger.debug(f"Starting {self.operation_name}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.start_time:
            duration = datetime.now().timestamp() - self.start_time
            self.logger.info(f"Completed {self.operation_name} in {duration:.3f}s")

            if exc_type:
                self.logger.error(f"Failed {self.operation_name}: {exc_val}")


def log_function_call(func):
    """Decorator to log function calls with timing."""

    def wrapper(*args, **kwargs):
        logger = get_analysis_logger()
        func_name = f"{func.__module__}.{func.__name__}"

        with PerformanceTimer(func_name, logger):
            try:
                result = func(*args, **kwargs)
                logger.debug(f"Successfully completed {func_name}")
                return result
            except Exception as e:
                logger.error(f"Error in {func_name}: {e}", exc_info=True)
                raise

    return wrapper


class LoggerMixin:
    """Mixin class to add logging capabilities to any class."""

    @property
    def logger(self) -> logging.Logger:
        """Get logger for this class."""
        return logging.getLogger(f"genoscope.{self.__class__.__name__}")

    def log_info(self, message: str, *args, **kwargs):
        """Log info message."""
        self.logger.info(message, *args, **kwargs)

    def log_warning(self, message: str, *args, **kwargs):
        """Log warning message."""
        self.logger.warning(message, *args, **kwargs)

    def log_error(self, message: str, *args, **kwargs):
        """Log error message."""
        self.logger.error(message, *args, **kwargs, exc_info=True)

    def log_debug(self, message: str, *args, **kwargs):
        """Log debug message."""
        self.logger.debug(message, *args, **kwargs)
