#!/bin/bash

echo "===================================="
echo "   Frontend Quick Fix & Start       "
echo "===================================="

cd /home/asd/.virtualenvs/BioForge_edit_branch/frontend

# Step 1: Clean everything
echo "🧹 Cleaning old dependencies..."
rm -rf node_modules package-lock.json yarn.lock
npm cache clean --force

# Step 2: Install dependencies
echo "📦 Installing dependencies (this will take 2-3 minutes)..."
npm install

# Step 3: Check if react-scripts is installed
if [ ! -f "node_modules/.bin/react-scripts" ]; then
    echo "⚠️  react-scripts not found, installing manually..."
    npm install react-scripts@5.0.1 --save-exact
fi

# Step 4: Try to start
echo "🚀 Starting frontend..."
npm start
