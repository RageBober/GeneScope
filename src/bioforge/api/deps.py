"""FastAPI dependencies."""

from typing import Annotated

from arq import ArqRedis, create_pool
from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from bioforge.config import Settings, get_settings
from bioforge.core.security import decode_token
from bioforge.core.exceptions import AuthenticationError
from bioforge.db import get_session
from bioforge.models import User


# ─────────────────────────────────────────────────────────────────────────────
# Type aliases for cleaner signatures
# ─────────────────────────────────────────────────────────────────────────────

SessionDep = Annotated[AsyncSession, Depends(get_session)]
SettingsDep = Annotated[Settings, Depends(get_settings)]


# ─────────────────────────────────────────────────────────────────────────────
# ARQ Task Queue
# ─────────────────────────────────────────────────────────────────────────────

_arq_pool: ArqRedis | None = None


async def get_arq_pool() -> ArqRedis:
    """Get ARQ Redis pool for enqueuing tasks."""
    global _arq_pool
    
    if _arq_pool is None:
        from bioforge.workers import get_redis_settings
        _arq_pool = await create_pool(get_redis_settings())
    
    return _arq_pool


async def close_arq_pool() -> None:
    """Close ARQ pool on shutdown."""
    global _arq_pool
    if _arq_pool is not None:
        await _arq_pool.close()
        _arq_pool = None


ArqPoolDep = Annotated[ArqRedis, Depends(get_arq_pool)]


# ─────────────────────────────────────────────────────────────────────────────
# Authentication
# ─────────────────────────────────────────────────────────────────────────────

security = HTTPBearer(auto_error=False)


async def get_current_user_optional(
    session: SessionDep,
    credentials: HTTPAuthorizationCredentials | None = Depends(security),
) -> User | None:
    """Get current user if authenticated, None otherwise."""
    if not credentials:
        return None
    
    try:
        payload = decode_token(credentials.credentials)
        user_id = payload.get("sub")
        if not user_id:
            return None
        
        result = await session.execute(
            select(User).where(User.id == int(user_id))
        )
        return result.scalar_one_or_none()
    except Exception:
        return None


async def get_current_user(
    user: User | None = Depends(get_current_user_optional),
) -> User:
    """Get current authenticated user. Raises 401 if not authenticated."""
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Not authenticated",
            headers={"WWW-Authenticate": "Bearer"},
        )
    return user


async def get_current_active_user(
    user: User = Depends(get_current_user),
) -> User:
    """Get current active user. Raises 403 if user is inactive."""
    if not user.is_active:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="User is inactive",
        )
    return user


async def get_current_superuser(
    user: User = Depends(get_current_active_user),
) -> User:
    """Get current superuser. Raises 403 if not superuser."""
    if not user.is_superuser:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Superuser access required",
        )
    return user


# Type aliases for user dependencies
CurrentUser = Annotated[User, Depends(get_current_user)]
CurrentActiveUser = Annotated[User, Depends(get_current_active_user)]
CurrentSuperuser = Annotated[User, Depends(get_current_superuser)]
OptionalUser = Annotated[User | None, Depends(get_current_user_optional)]
