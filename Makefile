.PHONY: install install-dev test lint format format-fix digest run nohup

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

# Install from source
install:
	pip install .

run:
	uv run main.py

nohup:
	bash -c 'nohup time uv run main.py & disown'

digest:
	gitingest -i "*.py, *.toml" . 
