#!/bin/bash

# GenoScope Frontend Start Script

echo "🚀 Starting GenoScope Frontend..."

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check if Node.js is installed
if ! command -v node &> /dev/null; then
    echo -e "${RED}❌ Node.js is not installed!${NC}"
    echo "Please install Node.js 16+ from: https://nodejs.org/"
    exit 1
fi

# Check Node version
NODE_VERSION=$(node -v | cut -d'v' -f2 | cut -d'.' -f1)
if [ "$NODE_VERSION" -lt 16 ]; then
    echo -e "${YELLOW}⚠️  Node.js version is less than 16. Please upgrade.${NC}"
fi

# Navigate to frontend directory
cd "$(dirname "$0")" || exit

# Install dependencies if node_modules doesn't exist
if [ ! -d "node_modules" ]; then
    echo -e "${YELLOW}📦 Installing dependencies...${NC}"
    npm install
    if [ $? -ne 0 ]; then
        echo -e "${RED}❌ Failed to install dependencies${NC}"
        exit 1
    fi
fi

# Create .env if it doesn't exist
if [ ! -f ".env" ]; then
    echo -e "${YELLOW}📝 Creating .env file...${NC}"
    cat > .env << EOL
REACT_APP_API_URL=http://localhost:8000/api
REACT_APP_WS_URL=ws://localhost:8000/ws
EOL
    echo -e "${GREEN}✅ .env file created${NC}"
fi

# Check if backend is running
echo -e "${YELLOW}🔍 Checking backend connection...${NC}"
if curl -s -o /dev/null -w "%{http_code}" http://localhost:8000/health | grep -q "200"; then
    echo -e "${GREEN}✅ Backend is running${NC}"
else
    echo -e "${YELLOW}⚠️  Backend is not running on port 8000${NC}"
    echo "Start the backend first with: cd .. && python src/main.py"
fi

# Start the development server
echo -e "${GREEN}🎉 Starting GenoScope Frontend...${NC}"
echo "-----------------------------------"
echo "📍 URL: http://localhost:3000"
echo "📍 API: http://localhost:8000"
echo "-----------------------------------"
echo -e "${YELLOW}Press Ctrl+C to stop${NC}"
echo ""

# Start React app
npm start
