#!/bin/bash

# GenoScope - Full Stack Startup Script

echo "🧬 Starting GenoScope Platform..."

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if a port is in use
port_in_use() {
    lsof -i:$1 >/dev/null 2>&1
}

# Function to start a service
start_service() {
    local name=$1
    local command=$2
    local port=$3
    
    echo -e "${BLUE}Starting $name...${NC}"
    
    if port_in_use $port; then
        echo -e "${YELLOW}Warning: Port $port is already in use. $name might already be running.${NC}"
    else
        eval "$command" &
        sleep 2
        if port_in_use $port; then
            echo -e "${GREEN}✓ $name started successfully on port $port${NC}"
        else
            echo -e "${RED}✗ Failed to start $name${NC}"
        fi
    fi
}

# Check prerequisites
echo -e "${BLUE}Checking prerequisites...${NC}"

if ! command_exists python3; then
    echo -e "${RED}Python 3 is not installed. Please install Python 3.9+${NC}"
    exit 1
fi

if ! command_exists node; then
    echo -e "${RED}Node.js is not installed. Please install Node.js 16+${NC}"
    exit 1
fi

if ! command_exists redis-server; then
    echo -e "${YELLOW}Redis is not installed. Some features may not work.${NC}"
fi

if ! command_exists psql; then
    echo -e "${YELLOW}PostgreSQL is not installed. Database features will not work.${NC}"
fi

echo -e "${GREEN}✓ Prerequisites check complete${NC}"

# Start Redis (if available)
if command_exists redis-server; then
    start_service "Redis" "redis-server" 6379
fi

# Start PostgreSQL (if available and not running)
if command_exists pg_ctl; then
    if ! pg_isready -q; then
        echo -e "${BLUE}Starting PostgreSQL...${NC}"
        pg_ctl start -D /usr/local/var/postgres -l /usr/local/var/postgres/server.log
        sleep 3
    else
        echo -e "${GREEN}✓ PostgreSQL is already running${NC}"
    fi
fi

# Navigate to project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Start Backend
echo -e "${BLUE}Starting Backend Services...${NC}"

# Activate Python virtual environment
if [ -f "venv/bin/activate" ]; then
    source venv/bin/activate
elif [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
else
    echo -e "${YELLOW}No virtual environment found. Creating one...${NC}"
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
fi

# Start FastAPI server
echo -e "${BLUE}Starting FastAPI server...${NC}"
uvicorn src.genoscope.api.main:app --reload --port 8000 --host 0.0.0.0 &
FASTAPI_PID=$!
sleep 3

# Start Celery worker
echo -e "${BLUE}Starting Celery worker...${NC}"
celery -A src.genoscope.api.tasks worker --loglevel=info &
CELERY_WORKER_PID=$!

# Start Celery beat
echo -e "${BLUE}Starting Celery beat...${NC}"
celery -A src.genoscope.api.tasks beat --loglevel=info &
CELERY_BEAT_PID=$!

# Start Frontend
echo -e "${BLUE}Starting Frontend...${NC}"
cd frontend

# Install dependencies if needed
if [ ! -d "node_modules" ]; then
    echo -e "${YELLOW}Installing frontend dependencies...${NC}"
    npm install
fi

# Start React development server
npm start &
FRONTEND_PID=$!

cd ..

# Wait for services to start
sleep 5

# Display status
echo ""
echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}🧬 GenoScope Platform Started Successfully!${NC}"
echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}"
echo ""
echo -e "${BLUE}Access Points:${NC}"
echo -e "  Frontend:    ${GREEN}http://localhost:3000${NC}"
echo -e "  Backend API: ${GREEN}http://localhost:8000${NC}"
echo -e "  API Docs:    ${GREEN}http://localhost:8000/docs${NC}"
echo -e "  GraphQL:     ${GREEN}http://localhost:8000/graphql${NC}"
echo ""
echo -e "${BLUE}Process IDs:${NC}"
echo -e "  FastAPI:     $FASTAPI_PID"
echo -e "  Celery Worker: $CELERY_WORKER_PID"
echo -e "  Celery Beat:   $CELERY_BEAT_PID"
echo -e "  Frontend:      $FRONTEND_PID"
echo ""
echo -e "${YELLOW}Press Ctrl+C to stop all services${NC}"

# Function to cleanup on exit
cleanup() {
    echo ""
    echo -e "${YELLOW}Shutting down GenoScope services...${NC}"
    
    # Kill all started processes
    kill $FASTAPI_PID 2>/dev/null
    kill $CELERY_WORKER_PID 2>/dev/null
    kill $CELERY_BEAT_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    
    echo -e "${GREEN}✓ All services stopped${NC}"
    exit 0
}

# Set up trap to cleanup on Ctrl+C
trap cleanup INT

# Keep script running
while true; do
    sleep 1
done
