#!/usr/bin/env python
"""
Простой скрипт для запуска тестов GenoScope
"""

import subprocess
import sys
from pathlib import Path

def run_command(cmd):
    """Запуск команды и вывод результата"""
    print(f"\n{'='*60}")
    print(f"Running: {cmd}")
    print('='*60)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr)
    return result.returncode

def main():
    base_dir = Path(__file__).parent
    
    print("\n" + "="*60)
    print("GenoScope Test Runner")
    print("="*60)
    
    # 1. Проверка структуры проекта
    print("\n1. Checking project structure...")
    if not (base_dir / "src").exists():
        print("❌ src/ directory not found")
    else:
        print("✅ src/ directory found")
    
    if not (base_dir / "tests").exists():
        print("❌ tests/ directory not found")
    else:
        print("✅ tests/ directory found")
    
    # 2. Запуск базовых тестов
    print("\n2. Running basic tests...")
    if (base_dir / "test_basic.py").exists():
        run_command(f"python {base_dir}/test_basic.py")
    
    # 3. Запуск рабочих тестов
    print("\n3. Running working tests...")
    if (base_dir / "test_working.py").exists():
        run_command(f"pytest {base_dir}/test_working.py -v")
    
    # 4. Проверка импортов
    print("\n4. Checking imports...")
    try:
        import fastapi
        print("✅ FastAPI installed")
    except ImportError:
        print("❌ FastAPI not installed")
    
    try:
        import pandas
        print("✅ Pandas installed")
    except ImportError:
        print("❌ Pandas not installed")
    
    try:
        import sqlalchemy
        print("✅ SQLAlchemy installed")
    except ImportError:
        print("❌ SQLAlchemy not installed")
    
    # 5. Проверка Docker
    print("\n5. Checking Docker...")
    docker_result = run_command("docker --version")
    if docker_result == 0:
        print("✅ Docker is available")
    else:
        print("❌ Docker is not available")
    
    # 6. Проверка наличия файлов
    print("\n6. Checking key files...")
    files_to_check = [
        "requirements.txt",
        "docker-compose.yml",
        "Dockerfile",
        "pyproject.toml",
        ".github/workflows/ci-cd.yml"
    ]
    
    for file in files_to_check:
        if (base_dir / file).exists():
            print(f"✅ {file} exists")
        else:
            print(f"❌ {file} not found")
    
    print("\n" + "="*60)
    print("Test run completed!")
    print("="*60)

if __name__ == "__main__":
    main()
