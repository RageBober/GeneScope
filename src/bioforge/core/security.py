"""Security utilities: file validation, authentication, JWT."""

import os
import re
from collections.abc import Iterable
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import aiofiles
from jose import JWTError, jwt
from passlib.context import CryptContext

from bioforge.config import settings
from bioforge.core.exceptions import (
    AuthenticationError,
    FileTooLargeError,
    InvalidFileTypeError,
)

# ─────────────────────────────────────────────────────────────────────────────
# Password Hashing
# ─────────────────────────────────────────────────────────────────────────────

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")


def hash_password(password: str) -> str:
    """Hash a password."""
    return pwd_context.hash(password)


def verify_password(plain_password: str, hashed_password: str) -> bool:
    """Verify a password against its hash."""
    return pwd_context.verify(plain_password, hashed_password)


# ─────────────────────────────────────────────────────────────────────────────
# JWT Tokens
# ─────────────────────────────────────────────────────────────────────────────

def create_access_token(
    data: dict[str, Any],
    expires_delta: timedelta | None = None,
) -> str:
    """Create a JWT access token."""
    to_encode = data.copy()
    expire = datetime.now(timezone.utc) + (
        expires_delta or timedelta(minutes=settings.security.access_token_expire_minutes)
    )
    to_encode.update({"exp": expire, "type": "access"})
    return jwt.encode(
        to_encode,
        settings.security.secret_key,
        algorithm=settings.security.algorithm,
    )


def create_refresh_token(data: dict[str, Any]) -> str:
    """Create a JWT refresh token."""
    to_encode = data.copy()
    expire = datetime.now(timezone.utc) + timedelta(
        days=settings.security.refresh_token_expire_days
    )
    to_encode.update({"exp": expire, "type": "refresh"})
    return jwt.encode(
        to_encode,
        settings.security.secret_key,
        algorithm=settings.security.algorithm,
    )


def decode_token(token: str) -> dict[str, Any]:
    """Decode and validate a JWT token."""
    try:
        payload = jwt.decode(
            token,
            settings.security.secret_key,
            algorithms=[settings.security.algorithm],
        )
        return payload
    except JWTError as e:
        raise AuthenticationError(f"Invalid token: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# File Security
# ─────────────────────────────────────────────────────────────────────────────

SAFE_CHARS = re.compile(r"[^A-Za-z0-9._-]+")


def secure_filename(name: str) -> str:
    """Sanitize filename for safe storage."""
    base = os.path.basename(name or "file")
    base = base.replace(" ", "_")
    base = SAFE_CHARS.sub("_", base)
    if not base or base in (".", ".."):
        base = "file"
    return base[:200]


def validate_file(
    path: Path,
    allowed_extensions: Iterable[str] | None = None,
    max_size: int | None = None,
) -> None:
    """
    Validate a saved file.
    
    Raises:
        FileNotFoundError: File doesn't exist
        FileTooLargeError: File exceeds size limit
        InvalidFileTypeError: File extension not allowed
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    stat = path.stat()
    
    if stat.st_size == 0:
        raise ValueError("Empty file")
    
    if max_size and stat.st_size > max_size:
        raise FileTooLargeError(stat.st_size, max_size)
    
    if allowed_extensions:
        # Handle double extensions like .fastq.gz
        ext = "".join(path.suffixes).lower()
        if ext not in allowed_extensions:
            # Try single extension
            ext = path.suffix.lower()
            if ext not in allowed_extensions:
                raise InvalidFileTypeError(ext, list(allowed_extensions))


def validate_saved_file(
    path: Path,
    allowed_exts: Iterable[str],
    max_size: int,
) -> tuple[bool, str]:
    """
    Legacy validation function.
    
    Returns (success, message) tuple for backward compatibility.
    """
    try:
        validate_file(path, allowed_exts, max_size)
        return True, "ok"
    except FileNotFoundError:
        return False, "File not found"
    except FileTooLargeError:
        return False, "File too large"
    except InvalidFileTypeError as e:
        return False, f"Extension {e.details.get('extension')} not allowed"
    except ValueError as e:
        return False, str(e)


async def async_save_upload(
    file,  # UploadFile
    dest: Path,
    max_size: int,
    chunk_size: int = 1024 * 1024,  # 1MB chunks
) -> int:
    """
    Save uploaded file asynchronously with size limit.
    
    Returns:
        Number of bytes written
        
    Raises:
        FileTooLargeError: If file exceeds max_size
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    bytes_written = 0
    
    async with aiofiles.open(dest, "wb") as out:
        while True:
            chunk = await file.read(chunk_size)
            if not chunk:
                break
            bytes_written += len(chunk)
            if bytes_written > max_size:
                await out.close()
                dest.unlink(missing_ok=True)
                raise FileTooLargeError(bytes_written, max_size)
            await out.write(chunk)
    
    return bytes_written
