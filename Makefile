.PHONY: install install-dev test lint format format-fix clean build install-local

# Install package in development mode
install:
	pip install -e .

# Install with development dependencies
install-dev:
	pip install -e ".[dev,viz]"

# Run tests
test:
	pytest

# Run linting
lint:
	ruff check foldfusion/ tests/ scripts/
	mypy foldfusion/

# Format code
format:
	ruff format foldfusion/ tests/ scripts/

# Format and fix linting issues
format-fix:
	ruff check --fix foldfusion/ tests/ scripts/
	ruff format foldfusion/ tests/ scripts/

# Clean build artifacts
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -type d -name __pycache__ -delete
	find . -name "*.pyc" -delete

# Build package
build:
	python -m build

# Install from source
install-local:
	pip install .

run:
	python3 main.py
