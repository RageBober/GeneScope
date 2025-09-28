@echo off
echo ===============================================
echo BioForge - Quick Fix and Start Script
echo ===============================================
echo.

echo [1] Applying automatic fixes...
python diagnostics\final_fix.py

echo.
echo [2] Installing dependencies...
pip install -r requirements.txt

echo.
echo [3] Testing imports...
python -c "from genoscope.main import main; print('OK: Main module')"
python -c "from genoscope.api.main import app; print('OK: API module')"
python -c "from genoscope.pipeline.qc import QCMetrics; print('OK: Pipeline module')"

echo.
echo [4] Starting API server...
echo.
echo ===============================================
echo Server will start at: http://localhost:8000
echo API Docs: http://localhost:8000/docs
echo Web UI: http://localhost:8000/ui
echo Press Ctrl+C to stop
echo ===============================================
echo.

cd src\genoscope\api
uvicorn main:app --reload --port 8000
