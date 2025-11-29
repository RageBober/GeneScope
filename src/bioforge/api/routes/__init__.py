"""API routes."""

from .health import router as health_router
from .auth import router as auth_router
from .projects import router as projects_router
from .samples import router as samples_router
from .jobs import router as jobs_router
from .variants import router as variants_router

__all__ = [
    "health_router",
    "auth_router",
    "projects_router",
    "samples_router",
    "jobs_router",
    "variants_router",
]
