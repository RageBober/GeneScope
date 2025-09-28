#!/usr/bin/env python3
"""Alternative script to run genoscope"""

import sys
import os

# Add src to Python path
project_root = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

try:
    from genoscope.main import main
    if __name__ == "__main__":
        main()
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure dependencies are installed: poetry install")
    sys.exit(1)
