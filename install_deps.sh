#!/bin/bash
# Install missing dependencies for GenoScope

echo "Installing missing dependencies..."

# Core dependencies
pip install aiofiles python-multipart httpx

# Optional but useful
pip install pytest-mock faker

echo "Dependencies installed!"
echo ""
echo "Now you can run tests with:"
echo "  pytest tests/ -v --no-cov"
echo "  or"
echo "  make test-quick"
