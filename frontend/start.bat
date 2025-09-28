@echo off
echo 🚀 Starting GenoScope Frontend...

:: Check if Node.js is installed
where node >nul 2>nul
if %errorlevel% neq 0 (
    echo ❌ Node.js is not installed!
    echo Please install Node.js 16+ from: https://nodejs.org/
    pause
    exit /b 1
)

:: Navigate to frontend directory
cd /d "%~dp0"

:: Install dependencies if node_modules doesn't exist
if not exist "node_modules" (
    echo 📦 Installing dependencies...
    call npm install
    if %errorlevel% neq 0 (
        echo ❌ Failed to install dependencies
        pause
        exit /b 1
    )
)

:: Create .env if it doesn't exist
if not exist ".env" (
    echo 📝 Creating .env file...
    (
        echo REACT_APP_API_URL=http://localhost:8000/api
        echo REACT_APP_WS_URL=ws://localhost:8000/ws
    ) > .env
    echo ✅ .env file created
)

:: Start the development server
echo.
echo 🎉 Starting GenoScope Frontend...
echo -----------------------------------
echo 📍 URL: http://localhost:3000
echo 📍 API: http://localhost:8000
echo -----------------------------------
echo Press Ctrl+C to stop
echo.

:: Start React app
call npm start
