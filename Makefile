# GenoScope Makefile
# Commands for project management

.PHONY: help install test test-all test-unit test-integration test-quick test-coverage clean run docker-up docker-down

help:
	@echo "GenoScope - Available commands:"
	@echo "  make install        - Install dependencies"
	@echo "  make test           - Run all tests"
	@echo "  make test-unit      - Run unit tests only"
	@echo "  make test-integration - Run integration tests"
	@echo "  make test-quick     - Run tests without coverage"
	@echo "  make test-coverage  - Run tests with coverage report"
	@echo "  make run            - Run the application"
	@echo "  make docker-up      - Start Docker containers"
	@echo "  make docker-down    - Stop Docker containers"
	@echo "  make clean          - Clean temporary files"

install:
	@echo "Installing dependencies..."
	pip install -r requirements.txt
	pip install -e .
	cd frontend && npm install

test:
	@echo "Running all tests..."
	pytest tests/ -v

test-all:
	@echo "Running all tests with coverage..."
	pytest tests/ -v --cov=src --cov-report=html --cov-report=term

test-unit:
	@echo "Running unit tests..."
	pytest tests/unit/ -v -m unit

test-integration:
	@echo "Running integration tests..."
	pytest tests/integration/ -v -m integration

test-quick:
	@echo "Running quick tests (no coverage)..."
	pytest tests/ -v --no-cov -x

test-coverage:
	@echo "Generating coverage report..."
	pytest tests/ --cov=src --cov-report=html --cov-report=term
	@echo "Coverage report generated in htmlcov/index.html"

test-framework:
	@echo "Testing framework setup..."
	pytest tests/test_framework.py -v

test-setup:
	@echo "Testing project setup..."
	pytest tests/unit/test_setup.py -v

run:
	@echo "Starting GenoScope..."
	python -m uvicorn src.genoscope.api.main:app --reload --host 0.0.0.0 --port 8000

run-frontend:
	@echo "Starting frontend..."
	cd frontend && npm start

docker-up:
	@echo "Starting Docker containers..."
	docker-compose up -d

docker-down:
	@echo "Stopping Docker containers..."
	docker-compose down

docker-build:
	@echo "Building Docker images..."
	docker-compose build

clean:
	@echo "Cleaning temporary files..."
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name ".coverage" -delete
	rm -rf htmlcov/
	rm -rf .pytest_cache/
	rm -rf test-results/
	rm -rf .coverage*

format:
	@echo "Formatting code..."
	black src/ tests/
	isort src/ tests/

lint:
	@echo "Running linters..."
	flake8 src/ tests/ --max-line-length=120
	mypy src/ --ignore-missing-imports

setup-dev:
	@echo "Setting up development environment..."
	pip install -e .
	pip install -r requirements.txt
	pip install pytest pytest-cov pytest-asyncio
	@echo "Development environment ready!"

check:
	@echo "Running all checks..."
	make lint
	make test-quick
	@echo "All checks passed!"

ci:
	@echo "Running CI pipeline..."
	make lint
	make test-all
	@echo "CI pipeline completed!"
