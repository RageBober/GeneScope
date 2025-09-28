#!/usr/bin/env python
"""
Test runner script for GenoScope
Run this from the project root directory
"""

import subprocess
import sys
from pathlib import Path

def run_tests():
    """Run all tests with proper configuration"""
    
    # Ensure we're in the project root
    project_root = Path(__file__).parent
    
    print("="*60)
    print("GenoScope Test Suite")
    print("="*60)
    print(f"Project root: {project_root}")
    print()
    
    # Different test commands
    commands = {
        "All tests": ["pytest", "tests/", "-v"],
        "Unit tests only": ["pytest", "tests/unit/", "-v", "-m", "unit"],
        "Integration tests": ["pytest", "tests/integration/", "-v", "-m", "integration"],
        "Quick tests (no coverage)": ["pytest", "tests/", "-v", "--no-cov"],
        "Test framework": ["pytest", "tests/test_framework.py", "-v"],
        "Setup tests": ["pytest", "tests/unit/test_setup.py", "-v"],
    }
    
    # Let user choose or run default
    if len(sys.argv) > 1:
        choice = sys.argv[1]
        if choice == "all":
            cmd = commands["All tests"]
        elif choice == "unit":
            cmd = commands["Unit tests only"]
        elif choice == "integration":
            cmd = commands["Integration tests"]
        elif choice == "quick":
            cmd = commands["Quick tests (no coverage)"]
        elif choice == "framework":
            cmd = commands["Test framework"]
        elif choice == "setup":
            cmd = commands["Setup tests"]
        else:
            print("Available options:")
            for key in commands.keys():
                print(f"  - {key.lower().replace(' ', '_')}")
            sys.exit(1)
    else:
        # Default: run all tests
        cmd = commands["All tests"]
    
    print(f"Running: {' '.join(cmd)}")
    print("-"*60)
    
    # Run the tests
    result = subprocess.run(cmd, cwd=project_root)
    
    return result.returncode

if __name__ == "__main__":
    sys.exit(run_tests())
