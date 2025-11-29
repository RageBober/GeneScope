"""BioForge entry point."""

import uvicorn

from bioforge.config import settings


def main():
    """Run the BioForge API server."""
    uvicorn.run(
        "bioforge.api.app:app",
        host=settings.host,
        port=settings.port,
        reload=settings.debug,
        log_level="debug" if settings.debug else "info",
    )


if __name__ == "__main__":
    main()
