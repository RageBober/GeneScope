#!/bin/bash

# GenoScope Platform Quick Start Script for Linux/WSL
# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Print colored messages
print_header() {
    echo -e "\n${BOLD}${BLUE}========================================${NC}"
    echo -e "${BOLD}${BLUE}   $1${NC}"
    echo -e "${BOLD}${BLUE}========================================${NC}\n"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if port is in use
port_in_use() {
    lsof -Pi :$1 -sTCP:LISTEN -t >/dev/null 2>&1
}

# Kill process on port
kill_port() {
    local port=$1
    if port_in_use $port; then
        print_warning "Port $port is in use. Attempting to free it..."
        kill -9 $(lsof -t -i:$port) 2>/dev/null
        sleep 1
        if port_in_use $port; then
            print_error "Could not free port $port"
            return 1
        else
            print_success "Port $port freed"
        fi
    fi
    return 0
}

# Cleanup function
cleanup() {
    print_info "Shutting down services..."
    
    # Kill backend
    if [ ! -z "$BACKEND_PID" ]; then
        kill $BACKEND_PID 2>/dev/null
        print_success "Backend stopped"
    fi
    
    # Kill frontend
    if [ ! -z "$FRONTEND_PID" ]; then
        kill $FRONTEND_PID 2>/dev/null
        print_success "Frontend stopped"
    fi
    
    # Kill any remaining processes on ports
    kill_port 8000 2>/dev/null
    kill_port 3000 2>/dev/null
    
    print_success "GenoScope Platform stopped"
    exit 0
}

# Trap Ctrl+C
trap cleanup INT

# Main script
print_header "GenoScope Platform - Quick Start"

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

print_info "Working directory: $SCRIPT_DIR"

# Check Python
print_info "Checking Python..."
if ! command_exists python3; then
    print_error "Python 3 is not installed"
    exit 1
fi

PYTHON_VERSION=$(python3 --version 2>&1 | grep -Po '(?<=Python )\d+\.\d+')
print_success "Python version: $PYTHON_VERSION"

# Check Node.js
print_info "Checking Node.js..."
if ! command_exists node; then
    print_error "Node.js is not installed. Please install Node.js 16+"
    print_info "Run: curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -"
    print_info "     sudo apt-get install -y nodejs"
    exit 1
fi

NODE_VERSION=$(node --version)
print_success "Node.js version: $NODE_VERSION"

# Setup Backend (located in src/genoscope)
print_header "Setting up Backend"
BACKEND_DIR="$SCRIPT_DIR/src/genoscope"

if [ ! -d "$BACKEND_DIR" ]; then
    print_error "Backend directory not found at $BACKEND_DIR"
    exit 1
fi

cd "$BACKEND_DIR"
print_success "Found backend at $BACKEND_DIR"

# Check for virtual environment in project root
VENV_DIR="$SCRIPT_DIR/.venv"
if [ ! -d "$VENV_DIR" ]; then
    VENV_DIR="$SCRIPT_DIR/venv"
    if [ ! -d "$VENV_DIR" ]; then
        print_info "Creating Python virtual environment..."
        python3 -m venv "$SCRIPT_DIR/.venv"
        VENV_DIR="$SCRIPT_DIR/.venv"
    fi
fi

print_info "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# Install backend dependencies
print_info "Installing backend dependencies..."

# Look for requirements.txt in different locations
if [ -f "$SCRIPT_DIR/requirements.txt" ]; then
    pip install -q -r "$SCRIPT_DIR/requirements.txt" 2>/dev/null || {
        print_warning "Some dependencies might be missing, but continuing..."
    }
elif [ -f "$BACKEND_DIR/requirements.txt" ]; then
    pip install -q -r "$BACKEND_DIR/requirements.txt" 2>/dev/null || {
        print_warning "Some dependencies might be missing, but continuing..."
    }
else
    print_warning "requirements.txt not found, assuming dependencies are installed"
fi

# Ensure essential packages are installed
pip install -q fastapi uvicorn sqlalchemy pydantic 2>/dev/null

print_success "Backend dependencies ready"

# Setup Frontend
print_header "Setting up Frontend"
FRONTEND_DIR="$SCRIPT_DIR/frontend"

if [ ! -d "$FRONTEND_DIR" ]; then
    print_error "Frontend directory not found at $FRONTEND_DIR"
    exit 1
fi

cd "$FRONTEND_DIR"
print_success "Found frontend at $FRONTEND_DIR"

if [ ! -f "package.json" ]; then
    print_error "package.json not found in frontend directory"
    exit 1
fi

if [ ! -d "node_modules" ]; then
    print_info "Installing frontend dependencies (this may take a few minutes)..."
    npm install --silent 2>/dev/null || {
        print_warning "Retrying with cache clean..."
        rm -rf node_modules package-lock.json
        npm cache clean --force
        npm install
    }
fi
print_success "Frontend dependencies ready"

# Start Backend
print_header "Starting Services"
cd "$BACKEND_DIR"

# Free port 8000 if needed
kill_port 8000

print_info "Starting backend server..."
source "$VENV_DIR/bin/activate"

# Start uvicorn from the genoscope directory
python3 -m uvicorn main:app --reload --host 0.0.0.0 --port 8000 > /tmp/backend.log 2>&1 &
BACKEND_PID=$!

# Wait for backend to start
sleep 3
if kill -0 $BACKEND_PID 2>/dev/null; then
    print_success "Backend started on http://localhost:8000"
    print_info "API docs: http://localhost:8000/docs"
else
    print_error "Failed to start backend"
    print_error "Check /tmp/backend.log for details:"
    tail -n 20 /tmp/backend.log
    exit 1
fi

# Start Frontend
cd "$FRONTEND_DIR"

# Free port 3000 if needed
kill_port 3000

print_info "Starting frontend server..."
PORT=3000 BROWSER=none npm start > /tmp/frontend.log 2>&1 &
FRONTEND_PID=$!

# Wait for frontend to start
print_info "Waiting for frontend to compile (this may take a minute)..."
sleep 10

if kill -0 $FRONTEND_PID 2>/dev/null; then
    print_success "Frontend started on http://localhost:3000"
else
    print_error "Failed to start frontend"
    print_error "Check /tmp/frontend.log for details:"
    tail -n 20 /tmp/frontend.log
    cleanup
    exit 1
fi

# Success message
print_header "GenoScope Platform is Running!"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}Frontend:${NC}    http://localhost:3000"
echo -e "${BOLD}Backend API:${NC} http://localhost:8000"
echo -e "${BOLD}API Docs:${NC}    http://localhost:8000/docs"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo
echo -e "${BOLD}Login Credentials:${NC}"
echo -e "  Email:    demo@genoscope.com"
echo -e "  Password: demo123"
echo
echo -e "${YELLOW}Press Ctrl+C to stop all services${NC}"
echo

# Keep running
while true; do
    if ! kill -0 $BACKEND_PID 2>/dev/null; then
        print_error "Backend stopped unexpectedly"
        break
    fi
    if ! kill -0 $FRONTEND_PID 2>/dev/null; then
        print_error "Frontend stopped unexpectedly"
        break
    fi
    sleep 2
done

cleanup
