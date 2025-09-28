#!/usr/bin/env python3
"""
GenoScope Platform - Quick Start Script
This script sets up and launches both backend and frontend services
"""

import os
import sys
import subprocess
import time
import signal
import psutil
from pathlib import Path

# Colors for terminal output
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'

def print_header(text):
    print(f"\n{Colors.BOLD}{Colors.HEADER}{'=' * 60}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER}{text.center(60)}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER}{'=' * 60}{Colors.END}\n")

def print_success(text):
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")

def print_error(text):
    print(f"{Colors.RED}✗ {text}{Colors.END}")

def print_info(text):
    print(f"{Colors.BLUE}ℹ {text}{Colors.END}")

def print_warning(text):
    print(f"{Colors.YELLOW}⚠ {text}{Colors.END}")

def check_command_exists(command):
    """Check if a command exists in the system"""
    return subprocess.call(f"which {command}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0

def check_port_available(port):
    """Check if a port is available"""
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) != 0

def kill_process_on_port(port):
    """Kill process running on specific port"""
    try:
        for proc in psutil.process_iter(['pid', 'name']):
            try:
                for conn in proc.connections():
                    if conn.laddr.port == port:
                        print_warning(f"Killing process {proc.info['name']} (PID: {proc.info['pid']}) on port {port}")
                        proc.kill()
                        time.sleep(1)
                        return True
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
    except:
        pass
    return False

def setup_backend():
    """Setup and check backend requirements"""
    print_header("Setting up Backend")
    
    # Backend is in src/genoscope
    backend_path = Path(__file__).parent / "src" / "genoscope"
    
    if not backend_path.exists():
        print_error(f"Backend directory not found at {backend_path}")
        return False
    
    os.chdir(backend_path)
    print_success(f"Found backend at {backend_path}")
    
    # Check Python version
    python_version = sys.version_info
    if python_version.major < 3 or python_version.minor < 8:
        print_error(f"Python 3.8+ required. Current version: {python_version.major}.{python_version.minor}")
        return False
    print_success(f"Python version: {python_version.major}.{python_version.minor}.{python_version.micro}")
    
    # Install backend dependencies
    print_info("Checking backend dependencies...")
    
    # Try to find requirements.txt in various locations
    requirements_locations = [
        backend_path / "requirements.txt",
        Path(__file__).parent / "requirements.txt",
        backend_path.parent.parent / "requirements.txt"
    ]
    
    requirements_file = None
    for location in requirements_locations:
        if location.exists():
            requirements_file = location
            print_success(f"Found requirements.txt at {requirements_file}")
            break
    
    if requirements_file:
        print_info("Installing backend dependencies...")
        result = subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-r", str(requirements_file)], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            print_warning("Some dependencies might not be installed, but continuing...")
        else:
            print_success("Backend dependencies installed")
    else:
        print_warning("requirements.txt not found, assuming dependencies are already installed")
    
    # Install additional required packages
    essential_packages = ["fastapi", "uvicorn", "sqlalchemy", "pydantic"]
    for package in essential_packages:
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", package], 
                      capture_output=True, text=True)
    
    return True

def setup_frontend():
    """Setup and check frontend requirements"""
    print_header("Setting up Frontend")
    
    frontend_path = Path(__file__).parent / "frontend"
    
    if not frontend_path.exists():
        print_error(f"Frontend directory not found at {frontend_path}")
        return False
    
    os.chdir(frontend_path)
    print_success(f"Found frontend at {frontend_path}")
    
    # Check Node.js
    if not check_command_exists("node"):
        print_error("Node.js is not installed. Please install Node.js 16+ first.")
        print_info("Visit: https://nodejs.org/")
        return False
    
    # Check Node version
    node_version = subprocess.run(["node", "--version"], capture_output=True, text=True)
    print_success(f"Node.js version: {node_version.stdout.strip()}")
    
    # Check npm
    if not check_command_exists("npm"):
        print_error("npm is not installed")
        return False
    
    npm_version = subprocess.run(["npm", "--version"], capture_output=True, text=True)
    print_success(f"npm version: {npm_version.stdout.strip()}")
    
    # Check if node_modules exists
    node_modules = frontend_path / "node_modules"
    package_json = frontend_path / "package.json"
    
    if not package_json.exists():
        print_error("package.json not found in frontend directory")
        return False
    
    if not node_modules.exists():
        print_info("Installing frontend dependencies (this may take a few minutes)...")
        result = subprocess.run(["npm", "install"], capture_output=True, text=True)
        if result.returncode != 0:
            print_error("Failed to install frontend dependencies")
            print(result.stderr)
            print_warning("Trying to clean cache and reinstall...")
            
            # Clean and retry
            subprocess.run(["rm", "-rf", "node_modules", "package-lock.json"])
            subprocess.run(["npm", "cache", "clean", "--force"])
            result = subprocess.run(["npm", "install"], capture_output=True, text=True)
            
            if result.returncode != 0:
                print_error("Failed to install dependencies after cleanup")
                return False
        
        print_success("Frontend dependencies installed")
    else:
        print_success("Frontend dependencies already installed")
    
    return True

def start_backend():
    """Start the backend server"""
    print_info("Starting backend server...")
    
    # Backend is in src/genoscope
    backend_path = Path(__file__).parent / "src" / "genoscope"
    os.chdir(backend_path)
    
    # Check if port 8000 is available
    if not check_port_available(8000):
        print_warning("Port 8000 is already in use")
        if kill_process_on_port(8000):
            print_success("Freed port 8000")
        else:
            print_error("Could not free port 8000, trying alternative port 8001")
            # Try alternative port
            process = subprocess.Popen(
                [sys.executable, "-m", "uvicorn", "main:app", "--reload", "--host", "0.0.0.0", "--port", "8001"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            time.sleep(3)
            if process.poll() is None:
                print_success("Backend server started on http://localhost:8001")
                print_info("API documentation available at http://localhost:8001/docs")
                return process, 8001
            else:
                print_error("Failed to start backend server")
                return None, None
    
    # Start FastAPI server
    process = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "main:app", "--reload", "--host", "0.0.0.0", "--port", "8000"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Wait for server to start
    time.sleep(3)
    
    if process.poll() is None:
        print_success("Backend server started on http://localhost:8000")
        print_info("API documentation available at http://localhost:8000/docs")
        return process, 8000
    else:
        print_error("Failed to start backend server")
        # Print error output for debugging
        stderr = process.stderr.read().decode()
        if stderr:
            print_error(f"Error: {stderr[:200]}")
        return None, None

def start_frontend():
    """Start the frontend server"""
    print_info("Starting frontend server...")
    
    frontend_path = Path(__file__).parent / "frontend"
    os.chdir(frontend_path)
    
    # Check if port 3000 is available
    port = 3000
    if not check_port_available(port):
        print_warning(f"Port {port} is already in use")
        if kill_process_on_port(port):
            print_success(f"Freed port {port}")
        else:
            # Try alternative port
            port = 3001
            if not check_port_available(port):
                print_error(f"Port {port} is also in use")
                return None, None
            print_info(f"Using alternative port {port}")
    
    # Start React development server
    env = os.environ.copy()
    env["PORT"] = str(port)
    env["BROWSER"] = "none"  # Don't auto-open browser
    
    process = subprocess.Popen(
        ["npm", "start"],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    
    # Wait for server to start
    print_info("Waiting for frontend to compile...")
    time.sleep(8)
    
    if process.poll() is None:
        print_success(f"Frontend server started on http://localhost:{port}")
        return process, port
    else:
        print_error("Failed to start frontend server")
        # Print error output for debugging
        stderr = process.stderr.read().decode()
        if stderr:
            print_error(f"Error: {stderr[:200]}")
        return None, None

def main():
    """Main function to orchestrate the startup"""
    print_header("GenoScope Platform - Quick Start")
    print_info("Starting GenoScope genomic analysis platform...")
    print_info(f"Working directory: {Path.cwd()}")
    
    # Change to project root
    project_root = Path(__file__).parent
    os.chdir(project_root)
    
    # Setup checks
    if not setup_backend():
        print_error("Backend setup failed")
        sys.exit(1)
    
    if not setup_frontend():
        print_error("Frontend setup failed")
        sys.exit(1)
    
    # Start services
    backend_process, backend_port = start_backend()
    if not backend_process:
        print_error("Could not start backend")
        sys.exit(1)
    
    frontend_process, frontend_port = start_frontend()
    if not frontend_process:
        print_error("Could not start frontend")
        if backend_process:
            backend_process.terminate()
        sys.exit(1)
    
    print_header("GenoScope Platform is Running!")
    print("")
    print("=" * 60)
    print_success(f"Frontend:    http://localhost:{frontend_port}")
    print_success(f"Backend API: http://localhost:{backend_port}")
    print_success(f"API Docs:    http://localhost:{backend_port}/docs")
    print("=" * 60)
    print()
    print_info("Login credentials:")
    print(f"  Email:    demo@genoscope.com")
    print(f"  Password: demo123")
    print()
    print_warning("Press Ctrl+C to stop all services")
    
    # Keep running and handle shutdown
    try:
        while True:
            time.sleep(1)
            
            # Check if processes are still running
            if backend_process.poll() is not None:
                print_error("Backend process stopped unexpectedly")
                stderr = backend_process.stderr.read().decode()
                if stderr:
                    print_error(f"Backend error: {stderr[:500]}")
                break
            if frontend_process.poll() is not None:
                print_error("Frontend process stopped unexpectedly")
                stderr = frontend_process.stderr.read().decode()
                if stderr:
                    print_error(f"Frontend error: {stderr[:500]}")
                break
                
    except KeyboardInterrupt:
        print("\n")
        print_info("Shutting down services...")
    finally:
        # Cleanup
        if backend_process and backend_process.poll() is None:
            backend_process.terminate()
            print_success("Backend stopped")
        
        if frontend_process and frontend_process.poll() is None:
            frontend_process.terminate()
            print_success("Frontend stopped")
        
        print_success("GenoScope Platform stopped")

if __name__ == "__main__":
    try:
        import psutil
    except ImportError:
        print_warning("Installing psutil for process management...")
        subprocess.run([sys.executable, "-m", "pip", "install", "psutil"])
        import psutil
    
    main()
