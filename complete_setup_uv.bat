@echo off
echo ===============================================
echo BioForge - Complete Setup with UV
echo ===============================================
echo.
echo This script will:
echo   1. Fix project issues
echo   2. Setup UV package manager
echo   3. Install all dependencies
echo   4. Start the API server
echo.
echo ===============================================
pause

echo.
echo [STEP 1] Fixing project issues...
echo ===============================================
python diagnostics\final_fix.py

echo.
echo [STEP 2] Installing UV package manager...
echo ===============================================
pip install uv

echo.
echo [STEP 3] Creating virtual environment with UV...
echo ===============================================
if exist .venv (
    echo Virtual environment already exists, skipping...
) else (
    uv venv
)

echo.
echo [STEP 4] Activating virtual environment...
echo ===============================================
call .venv\Scripts\activate

echo.
echo [STEP 5] Installing dependencies with UV (super fast!)...
echo ===============================================
uv pip install -r requirements.txt

echo.
echo [STEP 6] Installing development dependencies...
echo ===============================================
uv pip install pytest pytest-cov pytest-asyncio black ruff mypy ipython

echo.
echo [STEP 7] Installing project in editable mode...
echo ===============================================
uv pip install -e .

echo.
echo [STEP 8] Running health check...
echo ===============================================
python health_check.py

echo.
echo ===============================================
echo Setup Complete! 
echo ===============================================
echo.
echo Starting API server in 5 seconds...
echo.
echo Server will be available at:
echo   - API Docs: http://localhost:8000/docs
echo   - Web UI:   http://localhost:8000/ui
echo.
echo Press Ctrl+C to stop the server
echo ===============================================
timeout /t 5

cd src\genoscope\api
uvicorn main:app --reload --port 8000
