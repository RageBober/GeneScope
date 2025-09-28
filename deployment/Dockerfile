# GenoScope Production Dockerfile
# Multi-stage build for optimized production image

################################################################################
# Stage 1: Base Dependencies
################################################################################
FROM python:3.11-slim-bookworm AS base

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONHASHSEED=random \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    POETRY_VERSION=1.8.2 \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_CREATE=false \
    POETRY_CACHE_DIR=/tmp/poetry_cache \
    DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Essential build tools
    gcc \
    g++ \
    build-essential \
    make \
    # Libraries for bioinformatics packages
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    # System utilities
    curl \
    wget \
    git \
    # Clean up
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Create non-root user for security
RUN groupadd --gid 1000 genoscope \
    && useradd --uid 1000 --gid genoscope --shell /bin/bash --create-home genoscope

# Create application directories
WORKDIR /app
RUN mkdir -p /app/data /app/uploads /app/logs /app/temp \
    && chown -R genoscope:genoscope /app

################################################################################
# Stage 2: Poetry Dependencies
################################################################################
FROM base AS poetry

# Install Poetry
RUN pip install --no-cache-dir poetry==$POETRY_VERSION

# Copy dependency files
COPY pyproject.toml poetry.lock ./

# Configure Poetry and install dependencies
RUN poetry config virtualenvs.create false \
    && poetry install --only=main --extras="bio api" --no-dev --no-root \
    && rm -rf $POETRY_CACHE_DIR

################################################################################
# Stage 3: Application Build
################################################################################
FROM poetry AS builder

# Copy source code
COPY src/ ./src/
COPY README.md ./
COPY LICENSE ./

# Install the application
RUN poetry install --only-root

# Copy configuration and static assets
COPY config/ ./config/
COPY data/ ./data/

# Set correct permissions
RUN chown -R genoscope:genoscope /app

################################################################################
# Stage 4: Production Runtime
################################################################################
FROM python:3.11-slim-bookworm AS production

# Copy runtime environment variables from base
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONHASHSEED=random \
    PYTHONPATH=/app/src \
    DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Runtime libraries for bioinformatics packages
    libbz2-1.0 \
    liblzma5 \
    zlib1g \
    libcurl3-gnutls \
    libssl3 \
    # System utilities
    curl \
    # Health check utilities
    netcat-traditional \
    # Clean up
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Create non-root user
RUN groupadd --gid 1000 genoscope \
    && useradd --uid 1000 --gid genoscope --shell /bin/bash --create-home genoscope

# Set working directory
WORKDIR /app

# Copy Python dependencies from builder stage
COPY --from=builder /usr/local /usr/local

# Copy application files from builder stage
COPY --from=builder --chown=genoscope:genoscope /app /app

# Create additional directories for runtime
RUN mkdir -p /app/data/{datasets,reports,cache,backups} \
    /app/uploads \
    /app/logs \
    /app/temp \
    && chown -R genoscope:genoscope /app

# Switch to non-root user
USER genoscope

# Configure application settings
ENV GENOSCOPE_ENV=production \
    GENOSCOPE_HOST=0.0.0.0 \
    GENOSCOPE_PORT=8000 \
    LOG_LEVEL=INFO \
    UPLOAD_DIR=/app/uploads \
    DATA_DIR=/app/data \
    LOGS_DIR=/app/logs

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Expose port
EXPOSE 8000

# Default command
CMD ["uvicorn", "genoscope.api.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]

################################################################################
# Stage 5: Development Runtime (optional)
################################################################################
FROM production AS development

# Switch back to root for development tools installation
USER root

# Install development dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Development tools
    vim \
    nano \
    htop \
    iputils-ping \
    net-tools \
    # Development Python packages would be installed here if needed
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Install development Python dependencies
COPY --from=poetry $POETRY_CACHE_DIR $POETRY_CACHE_DIR
RUN pip install poetry==$POETRY_VERSION \
    && poetry config virtualenvs.create false

COPY pyproject.toml poetry.lock ./
RUN poetry install --extras="bio api cli all" --with dev

# Switch back to genoscope user
USER genoscope

# Override environment for development
ENV GENOSCOPE_ENV=development \
    LOG_LEVEL=DEBUG \
    GENOSCOPE_DEBUG=true

# Development command with auto-reload
CMD ["uvicorn", "genoscope.api.main:app", "--host", "0.0.0.0", "--port", "8000", "--reload", "--reload-dir", "/app/src"]

################################################################################
# Build Arguments and Labels
################################################################################

# Build arguments
ARG BUILD_DATE
ARG VCS_REF
ARG VERSION

# Labels for container metadata
LABEL maintainer="dzumenov@gmail.com" \
      org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="GenoScope" \
      org.label-schema.description="Advanced genomic data analysis platform" \
      org.label-schema.url="https://github.com/your-repo/genoscope" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/your-repo/genoscope" \
      org.label-schema.vendor="GenoScope Project" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0"