"""BioForge custom exceptions and error handlers."""

import logging
import traceback
from typing import Any

from fastapi import FastAPI, Request, status
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as StarletteHTTPException

logger = logging.getLogger("bioforge")


# ─────────────────────────────────────────────────────────────────────────────
# Custom Exceptions
# ─────────────────────────────────────────────────────────────────────────────

class BioForgeException(Exception):
    """Base exception for BioForge."""
    
    def __init__(
        self,
        message: str,
        code: str = "BIOFORGE_ERROR",
        status_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR,
        details: dict[str, Any] | None = None,
    ):
        self.message = message
        self.code = code
        self.status_code = status_code
        self.details = details or {}
        super().__init__(message)


class NotFoundError(BioForgeException):
    """Resource not found."""
    
    def __init__(self, resource: str, identifier: Any):
        super().__init__(
            message=f"{resource} not found: {identifier}",
            code="NOT_FOUND",
            status_code=status.HTTP_404_NOT_FOUND,
            details={"resource": resource, "identifier": str(identifier)},
        )


class ValidationError(BioForgeException):
    """Validation error."""
    
    def __init__(self, message: str, field: str | None = None):
        super().__init__(
            message=message,
            code="VALIDATION_ERROR",
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            details={"field": field} if field else {},
        )


class AuthenticationError(BioForgeException):
    """Authentication failed."""
    
    def __init__(self, message: str = "Authentication required"):
        super().__init__(
            message=message,
            code="AUTHENTICATION_ERROR",
            status_code=status.HTTP_401_UNAUTHORIZED,
        )


class AuthorizationError(BioForgeException):
    """Authorization failed."""
    
    def __init__(self, message: str = "Permission denied"):
        super().__init__(
            message=message,
            code="AUTHORIZATION_ERROR",
            status_code=status.HTTP_403_FORBIDDEN,
        )


class PipelineError(BioForgeException):
    """Pipeline execution error."""
    
    def __init__(self, message: str, step: str | None = None):
        super().__init__(
            message=message,
            code="PIPELINE_ERROR",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            details={"step": step} if step else {},
        )


class StorageError(BioForgeException):
    """Storage operation error."""
    
    def __init__(self, message: str, path: str | None = None):
        super().__init__(
            message=message,
            code="STORAGE_ERROR",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            details={"path": path} if path else {},
        )


class FileTooLargeError(BioForgeException):
    """File exceeds size limit."""
    
    def __init__(self, size: int, max_size: int):
        super().__init__(
            message=f"File too large: {size} bytes (max: {max_size})",
            code="FILE_TOO_LARGE",
            status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            details={"size": size, "max_size": max_size},
        )


class InvalidFileTypeError(BioForgeException):
    """Invalid file type."""
    
    def __init__(self, extension: str, allowed: list[str]):
        super().__init__(
            message=f"Invalid file type: {extension}",
            code="INVALID_FILE_TYPE",
            status_code=status.HTTP_415_UNSUPPORTED_MEDIA_TYPE,
            details={"extension": extension, "allowed": allowed},
        )


# ─────────────────────────────────────────────────────────────────────────────
# Exception Handlers
# ─────────────────────────────────────────────────────────────────────────────

def install_exception_handlers(app: FastAPI, include_trace: bool = False) -> None:
    """Install exception handlers on FastAPI app."""
    
    @app.exception_handler(BioForgeException)
    async def bioforge_exception_handler(request: Request, exc: BioForgeException):
        return JSONResponse(
            status_code=exc.status_code,
            content={
                "error": exc.code,
                "message": exc.message,
                "details": exc.details,
            },
        )
    
    @app.exception_handler(StarletteHTTPException)
    async def http_exception_handler(request: Request, exc: StarletteHTTPException):
        return JSONResponse(
            status_code=exc.status_code,
            content={
                "error": "HTTP_ERROR",
                "message": str(exc.detail),
            },
        )
    
    @app.exception_handler(RequestValidationError)
    async def validation_exception_handler(request: Request, exc: RequestValidationError):
        logger.warning("Validation error: %s", exc)
        return JSONResponse(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            content={
                "error": "VALIDATION_ERROR",
                "message": "Request validation failed",
                "details": exc.errors(),
            },
        )
    
    @app.exception_handler(Exception)
    async def unhandled_exception_handler(request: Request, exc: Exception):
        logger.exception("Unhandled error: %s", exc)
        content = {
            "error": "INTERNAL_ERROR",
            "message": "An unexpected error occurred",
        }
        if include_trace:
            content["trace"] = traceback.format_exc()
        return JSONResponse(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            content=content,
        )
