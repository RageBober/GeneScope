#!/bin/bash

# GenoScope - Stop All Services

echo "🧬 Stopping GenoScope Platform..."

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Kill processes by name
echo -e "${BLUE}Stopping FastAPI server...${NC}"
pkill -f "uvicorn.*genoscope"

echo -e "${BLUE}Stopping Celery workers...${NC}"
pkill -f "celery.*genoscope"

echo -e "${BLUE}Stopping React development server...${NC}"
pkill -f "react-scripts start"

echo -e "${BLUE}Stopping Node processes...${NC}"
pkill -f "node.*react"

# Optional: Stop Redis
if pgrep redis-server > /dev/null; then
    echo -e "${BLUE}Stopping Redis...${NC}"
    redis-cli shutdown
fi

echo ""
echo -e "${GREEN}✓ All GenoScope services stopped${NC}"
