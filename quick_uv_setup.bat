@echo off
echo ===============================================
echo BioForge - Quick UV Setup
echo ===============================================
echo.

echo [1] Installing UV...
pip install uv

echo.
echo [2] Creating virtual environment...
uv venv

echo.
echo [3] Activating environment...
call .venv\Scripts\activate

echo.
echo [4] Installing dependencies...
uv pip install -r requirements.txt

echo.
echo [5] Installing dev dependencies...
pip install pytest pytest-cov pytest-asyncio black ruff mypy

echo.
echo [6] Installing project in editable mode...
uv pip install -e .

echo.
echo ===============================================
echo UV Setup Complete!
echo ===============================================
echo.
echo To start the server:
echo   uvicorn genoscope.api.main:app --reload
echo.
echo To run tests:
echo   pytest tests/
echo.
echo UV Commands:
echo   uv pip list     - Show installed packages
echo   uv pip install  - Install a package
echo   uv pip freeze   - Show requirements
echo ===============================================
pause
