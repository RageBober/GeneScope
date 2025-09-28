#!/bin/bash

# Fix Frontend Installation Script
echo "===================================="
echo "   Fixing Frontend Dependencies     "
echo "===================================="

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Navigate to frontend directory
cd /home/asd/.virtualenvs/BioForge_edit_branch/frontend

echo -e "${YELLOW}Cleaning old installation...${NC}"
# Remove all old dependencies and cache
rm -rf node_modules package-lock.json yarn.lock

# Clean npm cache
npm cache clean --force

echo -e "${YELLOW}Installing dependencies...${NC}"
# Install dependencies with verbose output
npm install --verbose

# Check if react-scripts was installed
if [ -f "node_modules/.bin/react-scripts" ]; then
    echo -e "${GREEN}✓ react-scripts installed successfully${NC}"
else
    echo -e "${RED}✗ react-scripts not found, installing directly...${NC}"
    npm install react-scripts@5.0.1 --save
fi

# Verify all critical dependencies
echo -e "${YELLOW}Verifying critical packages...${NC}"

PACKAGES=("react" "react-dom" "react-scripts" "typescript" "react-router-dom" "@mui/material")

for package in "${PACKAGES[@]}"; do
    if [ -d "node_modules/$package" ]; then
        echo -e "${GREEN}✓ $package installed${NC}"
    else
        echo -e "${RED}✗ $package missing, installing...${NC}"
        npm install $package
    fi
done

echo -e "${GREEN}===================================="
echo -e "   Installation Complete!           "
echo -e "====================================${NC}"
echo ""
echo "Now try running: npm start"
