"""
Enhanced exception handling for GenoScope with custom exceptions and improved error responses.
"""

import logging
import traceback
from typing import Any

from fastapi import FastAPI
from fastapi import Request
from fastapi import status
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from pydantic import ValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException

logger = logging.getLogger("genoscope.exceptions")


# =============================================================================
# Custom Exception Classes
# =============================================================================


class GenoScopeException(Exception):
    """Base exception for all GenoScope-related errors."""

    def __init__(
        self,
        message: str,
        error_code: str = "GENOSCOPE_ERROR",
        details: dict | None = None,
    ):
        self.message = message
        self.error_code = error_code
        self.details = details or {}
        super().__init__(self.message)


class DataValidationError(GenoScopeException):
    """Raised when input data fails validation."""

    def __init__(
        self, message: str, field: str | None = None, details: dict | None = None
    ):
        super().__init__(message, "DATA_VALIDATION_ERROR", details)
        self.field = field


class FileProcessingError(GenoScopeException):
    """Raised when file processing fails."""

    def __init__(
        self, message: str, file_path: str | None = None, file_type: str | None = None
    ):
        details = {}
        if file_path:
            details["file_path"] = file_path
        if file_type:
            details["file_type"] = file_type
        super().__init__(message, "FILE_PROCESSING_ERROR", details)


class AnalysisError(GenoScopeException):
    """Raised when data analysis fails."""

    def __init__(
        self,
        message: str,
        analysis_type: str | None = None,
        details: dict | None = None,
    ):
        if analysis_type:
            details = details or {}
            details["analysis_type"] = analysis_type
        super().__init__(message, "ANALYSIS_ERROR", details)


class DatabaseError(GenoScopeException):
    """Raised when database operations fail."""

    def __init__(
        self, message: str, operation: str | None = None, table: str | None = None
    ):
        details = {}
        if operation:
            details["operation"] = operation
        if table:
            details["table"] = table
        super().__init__(message, "DATABASE_ERROR", details)


class ExternalServiceError(GenoScopeException):
    """Raised when external service calls fail."""

    def __init__(
        self,
        message: str,
        service_name: str | None = None,
        status_code: int | None = None,
    ):
        details = {}
        if service_name:
            details["service_name"] = service_name
        if status_code:
            details["status_code"] = status_code
        super().__init__(message, "EXTERNAL_SERVICE_ERROR", details)


class ConfigurationError(GenoScopeException):
    """Raised when configuration is invalid."""

    def __init__(self, message: str, config_key: str | None = None):
        details = {"config_key": config_key} if config_key else {}
        super().__init__(message, "CONFIGURATION_ERROR", details)


class ResourceNotFoundError(GenoScopeException):
    """Raised when a requested resource is not found."""

    def __init__(
        self,
        message: str,
        resource_type: str | None = None,
        resource_id: str | None = None,
    ):
        details = {}
        if resource_type:
            details["resource_type"] = resource_type
        if resource_id:
            details["resource_id"] = resource_id
        super().__init__(message, "RESOURCE_NOT_FOUND", details)


class InsufficientPermissionsError(GenoScopeException):
    """Raised when user lacks required permissions."""

    def __init__(self, message: str, required_permission: str | None = None):
        details = (
            {"required_permission": required_permission} if required_permission else {}
        )
        super().__init__(message, "INSUFFICIENT_PERMISSIONS", details)


# =============================================================================
# Exception Handlers
# =============================================================================


def create_error_response(
    error_code: str,
    message: str,
    status_code: int = 500,
    details: dict | None = None,
    include_trace: bool = False,
) -> JSONResponse:
    """Create standardized error response."""

    content = {
        "error": {
            "code": error_code,
            "message": message,
            "timestamp": __import__("datetime").datetime.utcnow().isoformat(),
        }
    }

    if details:
        content["error"]["details"] = details

    if include_trace:
        content["error"]["trace"] = traceback.format_exc()

    return JSONResponse(status_code=status_code, content=content)


def install_exception_handlers(app: FastAPI, include_trace: bool = False) -> None:
    """Install comprehensive exception handlers for the FastAPI app."""

    @app.exception_handler(GenoScopeException)
    async def genoscope_exception_handler(request: Request, exc: GenoScopeException):
        """Handle custom GenoScope exceptions."""
        logger.error(
            f"GenoScope error: {exc.message}",
            extra={
                "error_code": exc.error_code,
                "details": exc.details,
                "path": request.url.path,
                "method": request.method,
            },
        )

        # Map error codes to HTTP status codes
        status_map = {
            "DATA_VALIDATION_ERROR": status.HTTP_400_BAD_REQUEST,
            "FILE_PROCESSING_ERROR": status.HTTP_400_BAD_REQUEST,
            "ANALYSIS_ERROR": status.HTTP_422_UNPROCESSABLE_ENTITY,
            "DATABASE_ERROR": status.HTTP_500_INTERNAL_SERVER_ERROR,
            "EXTERNAL_SERVICE_ERROR": status.HTTP_503_SERVICE_UNAVAILABLE,
            "CONFIGURATION_ERROR": status.HTTP_500_INTERNAL_SERVER_ERROR,
            "RESOURCE_NOT_FOUND": status.HTTP_404_NOT_FOUND,
            "INSUFFICIENT_PERMISSIONS": status.HTTP_403_FORBIDDEN,
        }

        status_code = status_map.get(
            exc.error_code, status.HTTP_500_INTERNAL_SERVER_ERROR
        )

        return create_error_response(
            error_code=exc.error_code,
            message=exc.message,
            status_code=status_code,
            details=exc.details,
            include_trace=include_trace,
        )

    @app.exception_handler(StarletteHTTPException)
    async def http_exception_handler(request: Request, exc: StarletteHTTPException):
        """Handle HTTP exceptions."""
        logger.warning(
            f"HTTP error {exc.status_code}: {exc.detail}",
            extra={
                "path": request.url.path,
                "method": request.method,
                "status_code": exc.status_code,
            },
        )

        return create_error_response(
            error_code="HTTP_ERROR",
            message=exc.detail,
            status_code=exc.status_code,
            include_trace=include_trace,
        )

    @app.exception_handler(RequestValidationError)
    async def validation_exception_handler(
        request: Request, exc: RequestValidationError
    ):
        """Handle request validation errors."""
        logger.warning(
            f"Validation error: {exc.errors()}",
            extra={
                "path": request.url.path,
                "method": request.method,
                "errors": exc.errors(),
            },
        )

        # Format validation errors for better readability
        formatted_errors = []
        for error in exc.errors():
            field = " -> ".join(str(loc) for loc in error["loc"])
            formatted_errors.append(
                {"field": field, "message": error["msg"], "type": error["type"]}
            )

        return create_error_response(
            error_code="VALIDATION_ERROR",
            message="Request validation failed",
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            details={"errors": formatted_errors},
            include_trace=include_trace,
        )

    @app.exception_handler(ValidationError)
    async def pydantic_validation_handler(request: Request, exc: ValidationError):
        """Handle Pydantic validation errors."""
        logger.warning(f"Pydantic validation error: {exc.errors()}")

        return create_error_response(
            error_code="PYDANTIC_VALIDATION_ERROR",
            message="Data validation failed",
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            details={"errors": exc.errors()},
            include_trace=include_trace,
        )

    @app.exception_handler(Exception)
    async def unhandled_exception_handler(request: Request, exc: Exception):
        """Handle all unhandled exceptions."""
        logger.exception(
            "Unhandled exception occurred",
            extra={
                "path": request.url.path,
                "method": request.method,
                "exception_type": type(exc).__name__,
            },
        )

        return create_error_response(
            error_code="INTERNAL_SERVER_ERROR",
            message="An unexpected error occurred",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            details={"exception_type": type(exc).__name__},
            include_trace=include_trace,
        )


# =============================================================================
# Exception Context Managers and Decorators
# =============================================================================


class ExceptionContext:
    """Context manager for handling exceptions with logging and optional re-raising."""

    def __init__(
        self,
        operation_name: str,
        logger: logging.Logger | None = None,
        re_raise: bool = True,
        default_return: Any = None,
    ):
        self.operation_name = operation_name
        self.logger = logger or logging.getLogger("genoscope")
        self.re_raise = re_raise
        self.default_return = default_return

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            self.logger.error(
                f"Error in {self.operation_name}: {exc_val}", exc_info=True
            )

            if not self.re_raise:
                return True  # Suppress the exception
        return False


def handle_exceptions(
    operation_name: str | None = None,
    logger: logging.Logger | None = None,
    re_raise: bool = True,
    default_return: Any = None,
):
    """Decorator for handling exceptions in functions."""

    def decorator(func):
        def wrapper(*args, **kwargs):
            op_name = operation_name or f"{func.__module__}.{func.__name__}"
            func_logger = logger or logging.getLogger("genoscope")

            try:
                return func(*args, **kwargs)
            except Exception as e:
                func_logger.error(f"Error in {op_name}: {e}", exc_info=True)

                if re_raise:
                    raise
                return default_return

        return wrapper

    return decorator


def safe_execute(
    func,
    *args,
    operation_name: str | None = None,
    logger: logging.Logger | None = None,
    default_return: Any = None,
    **kwargs,
):
    """Safely execute a function with exception handling."""
    op_name = operation_name or f"{func.__name__}"
    exec_logger = logger or logging.getLogger("genoscope")

    try:
        return func(*args, **kwargs)
    except Exception as e:
        exec_logger.error(f"Error in {op_name}: {e}", exc_info=True)
        return default_return


# =============================================================================
# Validation Utilities
# =============================================================================


def validate_file_type(filename: str, allowed_extensions: set) -> None:
    """Validate file type against allowed extensions."""
    if not filename:
        raise DataValidationError("Filename is required")

    file_ext = filename.lower().split(".")[-1] if "." in filename else ""
    if f".{file_ext}" not in allowed_extensions:
        raise DataValidationError(
            f"File type '.{file_ext}' not allowed",
            field="filename",
            details={
                "allowed_extensions": list(allowed_extensions),
                "provided_extension": f".{file_ext}",
            },
        )


def validate_file_size(file_size: int, max_size: int) -> None:
    """Validate file size against maximum allowed size."""
    if file_size > max_size:
        raise DataValidationError(
            f"File size {file_size} bytes exceeds maximum allowed size",
            field="file_size",
            details={
                "file_size": file_size,
                "max_size": max_size,
                "size_mb": round(file_size / (1024 * 1024), 2),
                "max_mb": round(max_size / (1024 * 1024), 2),
            },
        )


def validate_numeric_range(
    value: float,
    min_value: float | None = None,
    max_value: float | None = None,
    field_name: str = "value",
) -> None:
    """Validate that a numeric value is within the specified range."""
    if min_value is not None and value < min_value:
        raise DataValidationError(
            f"{field_name} {value} is below minimum value {min_value}",
            field=field_name,
            details={"value": value, "min_value": min_value},
        )

    if max_value is not None and value > max_value:
        raise DataValidationError(
            f"{field_name} {value} is above maximum value {max_value}",
            field=field_name,
            details={"value": value, "max_value": max_value},
        )
