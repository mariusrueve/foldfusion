.PHONY: venv install install-dev test lint format format-fix digest run analysis nohup thesis thesis-clean thesis-lint docs docs-clean

VENV_DIR ?= .venv
PYTHON := $(VENV_DIR)/bin/python
PIP := $(VENV_DIR)/bin/pip

$(PYTHON):
	python3 -m venv $(VENV_DIR)
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install --upgrade wheel setuptools

venv: $(PYTHON)

# Install with development dependencies
install-dev: $(PYTHON)
	$(PYTHON) -m pip install -e ".[dev,viz]"

# Run tests
test: $(PYTHON)
	$(PYTHON) -m pytest

# Run linting
lint: $(PYTHON)
	$(PYTHON) -m ruff check foldfusion/ tests/ scripts/
	$(PYTHON) -m mypy foldfusion/

# Format code
format: $(PYTHON)
	$(PYTHON) -m ruff format foldfusion/ tests/ scripts/

# Format and fix linting issues
format-fix: $(PYTHON)
	$(PYTHON) -m ruff check --fix foldfusion/ tests/ scripts/
	$(PYTHON) -m ruff format foldfusion/ tests/ scripts/

# Install from source
install: $(PYTHON)
	$(PYTHON) -m pip install .

run: $(PYTHON)
	$(PYTHON) main.py

analysis:
	$(PYTHON) analysis/analysis.py

nohup: $(PYTHON)
	bash -c 'nohup $(PYTHON) main.py & disown'

digest:
	gitingest -i "*.py, *.toml" . 

# Build the LaTeX thesis with latexmk (biber configured in thesis/.latexmkrc)
thesis:
	cd thesis && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

# Build the operational documentation as a standalone PDF
docs:
	cd docs && latexmk -pdf -interaction=nonstopmode -halt-on-error documentation.tex

docs-clean:
	cd docs && latexmk -C

# Clean auxiliary LaTeX files via latexmk
thesis-clean:
	cd thesis && latexmk -C

# Lint LaTeX sources using chktex
thesis-lint:
	@if ! command -v chktex >/dev/null 2>&1; then \
		echo "error: chktex not found. Install it (e.g., via TeX Live) to lint LaTeX sources."; \
		exit 1; \
	fi
	cd thesis && find . -name '*.tex' -print0 | xargs -0 -r chktex -q
