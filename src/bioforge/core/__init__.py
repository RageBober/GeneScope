"""Core module - exceptions, security, logging, telemetry."""

from .exceptions import (
    BioForgeException,
    NotFoundError,
    ValidationError,
    AuthenticationError,
    AuthorizationError,
    PipelineError,
    StorageError,
    FileTooLargeError,
    InvalidFileTypeError,
    install_exception_handlers,
)
from .security import (
    hash_password,
    verify_password,
    create_access_token,
    create_refresh_token,
    decode_token,
    secure_filename,
    validate_file,
    validate_saved_file,
    async_save_upload,
)
from .logging import (
    setup_logging,
    get_logger,
    bind_context,
    clear_context,
)
from .telemetry import (
    setup_telemetry,
    instrument_fastapi,
    instrument_sqlalchemy,
    instrument_redis,
    get_tracer,
    traced,
)

__all__ = [
    # Exceptions
    "BioForgeException",
    "NotFoundError",
    "ValidationError",
    "AuthenticationError",
    "AuthorizationError",
    "PipelineError",
    "StorageError",
    "FileTooLargeError",
    "InvalidFileTypeError",
    "install_exception_handlers",
    # Security
    "hash_password",
    "verify_password",
    "create_access_token",
    "create_refresh_token",
    "decode_token",
    "secure_filename",
    "validate_file",
    "validate_saved_file",
    "async_save_upload",
    # Logging
    "setup_logging",
    "get_logger",
    "bind_context",
    "clear_context",
    # Telemetry
    "setup_telemetry",
    "instrument_fastapi",
    "instrument_sqlalchemy",
    "instrument_redis",
    "get_tracer",
    "traced",
]
