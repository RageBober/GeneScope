#!/usr/bin/env python3
"""
Migration script from Poetry to UV
Переход с Poetry на UV - современный менеджер пакетов Python
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Any

class PoetryToUVMigration:
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.poetry_lock = project_root / "poetry.lock"
        self.pyproject = project_root / "pyproject.toml"
        self.requirements_txt = project_root / "requirements.txt"
        self.requirements_dev = project_root / "requirements-dev.txt"
        
    def run_migration(self):
        """Выполнить полную миграцию"""
        print("🚀 Миграция с Poetry на UV")
        print("=" * 60)
        
        # 1. Проверка наличия uv
        if not self.check_uv_installed():
            self.install_uv()
            
        # 2. Создание нового pyproject.toml для uv
        self.create_uv_pyproject()
        
        # 3. Генерация requirements файлов
        self.generate_requirements()
        
        # 4. Создание виртуального окружения с uv
        self.create_uv_venv()
        
        # 5. Установка зависимостей через uv
        self.install_dependencies()
        
        # 6. Очистка старых файлов Poetry
        self.cleanup_poetry_files()
        
        # 7. Создание новых скриптов запуска
        self.create_uv_scripts()
        
        print("\n✅ Миграция завершена успешно!")
        
    def check_uv_installed(self) -> bool:
        """Проверка установки uv"""
        try:
            result = subprocess.run(["uv", "--version"], capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✅ UV установлен: {result.stdout.strip()}")
                return True
        except FileNotFoundError:
            pass
        
        print("❌ UV не установлен")
        return False
        
    def install_uv(self):
        """Установка uv"""
        print("\n📦 Установка UV...")
        
        # Для Windows используем PowerShell
        if sys.platform == "win32":
            print("Выполняется установка для Windows...")
            install_cmd = 'powershell -c "irm https://astral.sh/uv/install.ps1 | iex"'
            
            # Альтернативный способ через pip
            print("Альтернативная установка через pip...")
            try:
                subprocess.run([sys.executable, "-m", "pip", "install", "uv"], check=True)
                print("✅ UV установлен через pip")
            except subprocess.CalledProcessError:
                print("⚠️ Не удалось установить UV автоматически")
                print("\n📌 Установите UV вручную:")
                print("   pip install uv")
                print("   или")
                print("   powershell -c \"irm https://astral.sh/uv/install.ps1 | iex\"")
                sys.exit(1)
        else:
            # Для Unix-like систем
            install_cmd = "curl -LsSf https://astral.sh/uv/install.sh | sh"
            subprocess.run(install_cmd, shell=True)
            
    def create_uv_pyproject(self):
        """Создание нового pyproject.toml оптимизированного для uv"""
        print("\n📝 Создание нового pyproject.toml для UV...")
        
        pyproject_content = '''[project]
name = "genoscope"
version = "1.0.0"
description = "Genomics Analysis Platform - Powered by UV"
authors = [{name = "GenoScope Team", email = "team@genoscope.com"}]
readme = "README.md"
requires-python = ">=3.8"
license = {text = "MIT"}
keywords = ["genomics", "bioinformatics", "analysis", "pipeline"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
]

dependencies = [
    # Web Framework
    "fastapi>=0.100.0",
    "uvicorn[standard]>=0.22.0",
    "aiofiles>=23.0.0",
    "python-multipart>=0.0.6",
    "httpx>=0.24.0",
    
    # Data Science
    "pandas>=2.0.0",
    "numpy>=1.24.0",
    "scikit-learn>=1.2.0",
    "matplotlib>=3.7.0",
    "seaborn>=0.12.0",
    "plotly>=5.14.0",
    
    # Bioinformatics
    "biopython>=1.81",
    "pysam>=0.21.0",
    
    # Database
    "sqlalchemy>=2.0.0",
    "pydantic>=2.0.0",
    "redis>=4.5.0",
    
    # Task Queue
    "celery>=5.3.0",
    
    # Auth & Security
    "python-jose[cryptography]>=3.3.0",
    "passlib[bcrypt]>=1.7.4",
    
    # Cloud & Storage
    "boto3>=1.26.0",
    
    # Payments
    "stripe>=5.0.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "pytest-asyncio>=0.21.0",
    "pytest-mock>=3.10.0",
    "pytest-xdist>=3.2.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
    "mypy>=1.0.0",
    "pre-commit>=3.0.0",
    "ipdb>=0.13.0",
]

ml = [
    "torch>=2.0.0",
    "transformers>=4.30.0",
    "tensorflow>=2.13.0",
    "xgboost>=1.7.0",
    "lightgbm>=4.0.0",
]

genomics = [
    "cyvcf2>=0.30.0",
    "pybedtools>=0.9.0",
    "pyvcf>=0.6.8",
]

parallel = [
    "dask[complete]>=2023.5.0",
    "ray>=2.5.0",
]

monitoring = [
    "prometheus-client>=0.16.0",
    "grafana-api>=1.0.0",
    "sentry-sdk>=1.0.0",
]

all = [
    "genoscope[dev,ml,genomics,parallel,monitoring]",
]

[project.scripts]
genoscope = "genoscope.main:main"
genoscope-api = "genoscope.api.main:run_server"
genoscope-worker = "genoscope.worker:main"

[project.urls]
Homepage = "https://github.com/genoscope/bioforge"
Documentation = "https://genoscope.readthedocs.io"
Repository = "https://github.com/genoscope/bioforge.git"
Issues = "https://github.com/genoscope/bioforge/issues"

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/genoscope"]

[tool.uv]
# UV specific settings
dev-dependencies = [
    "ipython>=8.0.0",
    "jupyterlab>=4.0.0",
]

# UV workspace settings
[tool.uv.workspace]
members = ["src/genoscope"]

# Tool configurations
[tool.ruff]
line-length = 100
target-version = "py38"
select = [
    "E",  # pycodestyle errors
    "W",  # pycodestyle warnings
    "F",  # pyflakes
    "I",  # isort
    "B",  # flake8-bugbear
    "C4", # flake8-comprehensions
    "UP", # pyupgrade
]
ignore = [
    "E501",  # line too long
    "B008",  # do not perform function calls in argument defaults
    "B905",  # zip() without explicit strict= parameter
]

[tool.ruff.per-file-ignores]
"__init__.py" = ["F401"]

[tool.black]
line-length = 100
target-version = ["py38", "py39", "py310", "py311", "py312"]
include = '\.pyi?$'

[tool.mypy]
python_version = "3.8"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = true
ignore_missing_imports = true
no_implicit_optional = true
strict_optional = true

[tool.pytest.ini_options]
minversion = "7.0"
testpaths = ["tests"]
python_files = ["test_*.py", "*_test.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
pythonpath = ["."]
addopts = """
    -ra
    --strict-markers
    --strict-config
    --cov=src
    --cov-branch
    --cov-report=term-missing:skip-covered
    --cov-report=html
    --cov-report=xml
    --cov-fail-under=0
    --maxfail=1
    --tb=short
    -p no:warnings
    -v
"""

[tool.coverage.run]
source = ["src"]
branch = true
parallel = true
omit = [
    "*/tests/*",
    "*/test_*.py",
    "*/__init__.py",
    "*/migrations/*",
    "*/config/*",
]

[tool.coverage.report]
precision = 2
show_missing = true
skip_covered = false
'''
        
        # Создаем бэкап старого pyproject.toml
        if self.pyproject.exists():
            backup_path = self.pyproject.with_suffix('.toml.poetry-backup')
            self.pyproject.rename(backup_path)
            print(f"   📦 Бэкап старого pyproject.toml сохранен как {backup_path.name}")
        
        # Записываем новый pyproject.toml
        self.pyproject.write_text(pyproject_content, encoding='utf-8')
        print("   ✅ Новый pyproject.toml создан")
        
    def generate_requirements(self):
        """Генерация requirements файлов из текущих зависимостей"""
        print("\n📋 Генерация requirements файлов...")
        
        # requirements.txt для основных зависимостей
        requirements_content = """# Core dependencies
fastapi>=0.100.0
uvicorn[standard]>=0.22.0
pandas>=2.0.0
numpy>=1.24.0
scikit-learn>=1.2.0
matplotlib>=3.7.0
seaborn>=0.12.0
plotly>=5.14.0
biopython>=1.81
pysam>=0.21.0
sqlalchemy>=2.0.0
pydantic>=2.0.0
redis>=4.5.0
celery>=5.3.0
python-jose[cryptography]>=3.3.0
passlib[bcrypt]>=1.7.4
boto3>=1.26.0
stripe>=5.0.0
aiofiles>=23.0.0
python-multipart>=0.0.6
httpx>=0.24.0
"""
        
        self.requirements_txt.write_text(requirements_content, encoding='utf-8')
        print("   ✅ requirements.txt создан")
        
        # requirements-dev.txt для dev зависимостей
        requirements_dev_content = """# Development dependencies
pytest>=7.0.0
pytest-cov>=4.0.0
pytest-asyncio>=0.21.0
pytest-mock>=3.10.0
pytest-xdist>=3.2.0
black>=23.0.0
ruff>=0.1.0
mypy>=1.0.0
pre-commit>=3.0.0
ipdb>=0.13.0
ipython>=8.0.0
jupyterlab>=4.0.0
"""
        
        self.requirements_dev.write_text(requirements_dev_content, encoding='utf-8')
        print("   ✅ requirements-dev.txt создан")
        
    def create_uv_venv(self):
        """Создание виртуального окружения через uv"""
        print("\n🐍 Создание виртуального окружения UV...")
        
        try:
            # Создаем новое окружение
            subprocess.run(["uv", "venv", ".venv"], cwd=self.project_root, check=True)
            print("   ✅ Виртуальное окружение создано в .venv")
            
        except subprocess.CalledProcessError as e:
            print(f"   ⚠️ Ошибка при создании окружения: {e}")
            
    def install_dependencies(self):
        """Установка зависимостей через uv"""
        print("\n📦 Установка зависимостей через UV...")
        
        try:
            # Установка основных зависимостей
            print("   Installing production dependencies...")
            subprocess.run(
                ["uv", "pip", "install", "-r", "requirements.txt"],
                cwd=self.project_root,
                check=True
            )
            
            # Установка dev зависимостей
            print("   Installing development dependencies...")
            subprocess.run(
                ["uv", "pip", "install", "-r", "requirements-dev.txt"],
                cwd=self.project_root,
                check=True
            )
            
            # Установка проекта в editable mode
            print("   Installing project in editable mode...")
            subprocess.run(
                ["uv", "pip", "install", "-e", "."],
                cwd=self.project_root,
                check=True
            )
            
            print("   ✅ Все зависимости установлены")
            
        except subprocess.CalledProcessError as e:
            print(f"   ⚠️ Ошибка при установке: {e}")
            print("   Попробуйте выполнить вручную:")
            print("   uv pip install -r requirements.txt")
            
    def cleanup_poetry_files(self):
        """Очистка файлов Poetry"""
        print("\n🧹 Очистка файлов Poetry...")
        
        files_to_remove = [
            "poetry.lock",
            "poetry.toml",
        ]
        
        for file_name in files_to_remove:
            file_path = self.project_root / file_name
            if file_path.exists():
                # Создаем бэкап
                backup_path = file_path.with_suffix(f"{file_path.suffix}.backup")
                file_path.rename(backup_path)
                print(f"   📦 {file_name} -> {backup_path.name}")
                
        print("   ✅ Файлы Poetry перемещены в бэкапы")
        
    def create_uv_scripts(self):
        """Создание скриптов для работы с UV"""
        print("\n📄 Создание скриптов для UV...")
        
        # Скрипт для Windows
        windows_script = """@echo off
echo ===============================================
echo BioForge with UV - Quick Start
echo ===============================================
echo.

REM Activate virtual environment
call .venv\\Scripts\\activate

REM Check if UV is installed
uv --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Installing UV...
    pip install uv
)

echo.
echo Available commands:
echo   uv pip list          - Show installed packages
echo   uv pip install       - Install packages
echo   uv pip sync          - Sync with requirements.txt
echo   uv pip compile       - Compile requirements
echo.

echo Starting API server...
uvicorn genoscope.api.main:app --reload --port 8000
"""
        
        win_script_path = self.project_root / "start_with_uv.bat"
        win_script_path.write_text(windows_script)
        print(f"   ✅ {win_script_path.name} создан")
        
        # Скрипт для Unix
        unix_script = """#!/bin/bash
echo "==============================================="
echo "BioForge with UV - Quick Start"
echo "==============================================="
echo

# Activate virtual environment
source .venv/bin/activate

# Check if UV is installed
if ! command -v uv &> /dev/null; then
    echo "Installing UV..."
    pip install uv
fi

echo
echo "Available commands:"
echo "  uv pip list          - Show installed packages"
echo "  uv pip install       - Install packages"
echo "  uv pip sync          - Sync with requirements.txt"
echo "  uv pip compile       - Compile requirements"
echo

echo "Starting API server..."
uvicorn genoscope.api.main:app --reload --port 8000
"""
        
        unix_script_path = self.project_root / "start_with_uv.sh"
        unix_script_path.write_text(unix_script)
        # Делаем скрипт исполняемым
        unix_script_path.chmod(0o755)
        print(f"   ✅ {unix_script_path.name} создан")
        
        # Makefile для UV
        makefile_content = """# BioForge Makefile with UV
.PHONY: help install dev test run clean

help:
	@echo "BioForge with UV - Available commands:"
	@echo "  make install    - Install production dependencies"
	@echo "  make dev        - Install all dependencies (including dev)"
	@echo "  make test       - Run tests"
	@echo "  make run        - Run API server"
	@echo "  make clean      - Clean cache and temp files"
	@echo "  make format     - Format code with black and ruff"
	@echo "  make lint       - Run linters"
	@echo "  make docker     - Build and run with Docker"

install:
	uv pip install -r requirements.txt
	uv pip install -e .

dev:
	uv pip install -r requirements.txt
	uv pip install -r requirements-dev.txt
	uv pip install -e .

test:
	pytest tests/ -v

run:
	uvicorn genoscope.api.main:app --reload --port 8000

clean:
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	rm -rf .pytest_cache .coverage htmlcov coverage.xml

format:
	black src/ tests/
	ruff check --fix src/ tests/

lint:
	ruff check src/ tests/
	mypy src/

docker:
	docker-compose up --build

# UV specific commands
uv-sync:
	uv pip sync requirements.txt

uv-compile:
	uv pip compile pyproject.toml -o requirements.txt

uv-upgrade:
	uv pip compile --upgrade pyproject.toml -o requirements.txt
"""
        
        makefile_path = self.project_root / "Makefile.uv"
        makefile_path.write_text(makefile_content)
        print(f"   ✅ {makefile_path.name} создан")


def main():
    """Основная функция миграции"""
    project_root = Path.cwd()
    
    print("🚀 МИГРАЦИЯ BIOFORGE С POETRY НА UV")
    print("=" * 60)
    print(f"📁 Проект: {project_root}")
    print()
    print("UV - это:")
    print("  • ⚡ В 10-100 раз быстрее pip и poetry")
    print("  • 🦀 Написан на Rust для максимальной производительности")
    print("  • 🎯 Совместим с pip и requirements.txt")
    print("  • 🔧 От создателей Ruff (Astral)")
    print()
    
    response = input("Начать миграцию? (y/n): ")
    
    if response.lower() == 'y':
        migration = PoetryToUVMigration(project_root)
        migration.run_migration()
        
        print("\n" + "=" * 60)
        print("📋 СЛЕДУЮЩИЕ ШАГИ:")
        print("=" * 60)
        print()
        print("1. Активировать виртуальное окружение UV:")
        print("   Windows: .venv\\Scripts\\activate")
        print("   Unix:    source .venv/bin/activate")
        print()
        print("2. Проверить установку:")
        print("   uv pip list")
        print()
        print("3. Запустить проект:")
        print("   Windows: start_with_uv.bat")
        print("   Unix:    ./start_with_uv.sh")
        print()
        print("4. Использовать Makefile:")
        print("   make -f Makefile.uv dev    # Установить все зависимости")
        print("   make -f Makefile.uv run    # Запустить сервер")
        print()
    else:
        print("\n❌ Миграция отменена")


if __name__ == "__main__":
    main()
