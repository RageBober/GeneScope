#!/bin/bash

echo "🔧 Creating new React app from scratch..."

# Navigate to parent directory
cd /home/asd/.virtualenvs/BioForge_edit_branch/

# Remove old frontend if exists
if [ -d "frontend" ]; then
    echo "Removing old frontend directory..."
    rm -rf frontend
fi

# Create new React app with TypeScript
echo "Creating new React app with TypeScript..."
npx create-react-app frontend --template typescript

cd frontend

# Install all required dependencies
echo "Installing Material-UI and other dependencies..."
npm install @mui/material@^5.14.18 @emotion/react@^11.11.1 @emotion/styled@^11.11.0
npm install @mui/icons-material@^5.14.18 @mui/lab@^5.0.0-alpha.153
npm install react-router-dom@^6.20.0
npm install axios@^1.6.2
npm install lucide-react@^0.292.0
npm install notistack@^3.0.1
npm install react-dropzone@^14.2.3
npm install recharts@^2.9.3
npm install date-fns@^2.30.0

echo "✅ New frontend created successfully!"
echo "Now copy the src files from the backup"
