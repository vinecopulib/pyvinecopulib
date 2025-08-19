# Makefile for pyvinecopulib development
.PHONY: help install install-dev install-docs build clean test test-fast test-cov test-examples coverage coverage-open lint format type-check security docs docs-serve pre-commit update-deps check-all release
.DEFAULT_GOAL := help

# Python and package manager commands
# Prefer conda/mamba if available, fallback to pip
CONDA := $(shell command -v mamba 2> /dev/null || command -v conda 2> /dev/null)
PYTHON := python
PIP := pip

# Check if we're in a conda environment
ifdef CONDA_DEFAULT_ENV
	PKG_INSTALL := $(CONDA) install -c conda-forge
	PIP_INSTALL := $(PIP) install
else
	PKG_INSTALL := $(PIP) install
	PIP_INSTALL := $(PIP) install
endif

# Project directories
SRC_DIR := src
TEST_DIR := tests
DOCS_DIR := docs
EXAMPLES_DIR := examples

help: ## Show this help message
	@echo "pyvinecopulib development commands:"
	@echo
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install the package in development mode
	$(PIP_INSTALL) -e .

install-dev: ## Install development dependencies
	$(PIP_INSTALL) -e ".[dev]"

install-docs: ## Install documentation dependencies
	$(PIP_INSTALL) -e ".[doc]"

install-examples: ## Install examples dependencies
	$(PIP_INSTALL) -e ".[examples]"

install-all: ## Install all dependencies (dev, docs, examples)
	$(PIP_INSTALL) -e ".[dev,doc,examples]"

clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete

test: ## Run all tests (with coverage if installed)
	$(PYTHON) -m pytest $(TEST_DIR) -v

test-fast: ## Run tests without coverage
	$(PYTHON) -m pytest $(TEST_DIR) -x --no-cov

test-examples: ## Run example notebooks
	$(PYTHON) -m pytest --nbmake $(EXAMPLES_DIR)

coverage: ## Generate coverage report only (run tests first)
	$(PYTHON) -m coverage report --show-missing
	$(PYTHON) -m coverage html

coverage-open: ## Open coverage report in browser
	$(PYTHON) -c "import webbrowser; webbrowser.open('htmlcov/index.html')"

lint: ## Run linting with ruff
	$(PYTHON) -m ruff check $(SRC_DIR) $(TEST_DIR) --fix

format: ## Format code with ruff
	$(PYTHON) -m ruff format $(SRC_DIR) $(TEST_DIR)

type-check: ## Run type checking with mypy
	$(PYTHON) -m mypy

docs: ## Build documentation
	cd $(DOCS_DIR) && $(PYTHON) -m sphinx -b html . _build/html

docs-serve: ## Serve documentation locally
	cd $(DOCS_DIR) && $(PYTHON) serve_sphinx.py

docs-clean: ## Clean documentation build
	cd $(DOCS_DIR) && rm -rf _build/

pre-commit-install: ## Install pre-commit hooks
	$(PIP_INSTALL) pre-commit
	pre-commit install

pre-commit: ## Run pre-commit on all files
	pre-commit run --all-files

update-deps: ## Update dependencies
	$(PYTHON) scripts/generate_requirements.py --format txt
	$(PYTHON) scripts/generate_requirements.py --format yml

clear-cache: ## Clear Python cache files
	zsh scripts/clear_cache.sh

check-all: lint type-check security test-cov ## Run all quality checks

stubs: ## Generate type stubs using custom script
	$(PYTHON) scripts/generate_metadata.py --env $(CONDA_DEFAULT_ENV) --no-docstrings --no-examples

docstrings: ## Generate C++ docstrings
	$(PYTHON) scripts/generate_metadata.py --env $(CONDA_DEFAULT_ENV) --no-stubs --no-examples

metadata: ## Generate all metadata (docstrings, stubs, examples)
	$(PYTHON) scripts/generate_metadata.py --env $(CONDA_DEFAULT_ENV)

examples: ## Process and execute example notebooks
	$(PYTHON) scripts/generate_metadata.py --env $(CONDA_DEFAULT_ENV) --no-docstrings --no-stubs

# Development workflow commands
dev-setup: install-all pre-commit-install ## Complete development setup
	@echo "Development environment setup complete!"
	@echo "Run 'make help' to see available commands."

quick-check: lint type-check test-fast ## Quick development check (fast)

# Release workflow
release-check: clean check-all test-examples docs ## Pre-release checks
	@echo "All release checks passed!"

# Environment management
env-conda: ## Create conda environment
	$(PYTHON) scripts/generate_requirements.py --format yml
	$(CONDA) env create -f environment.yml

env-update: ## Update conda environment
	$(PYTHON) scripts/generate_requirements.py --format yml
	$(CONDA) env update -f environment.yml

env-activate: ## Show command to activate conda environment
	@echo "Run: conda activate pyvinecopulib"

# Git helpers
git-clean: ## Clean git working directory
	git clean -fdx

status: ## Show project status
	@echo "=== Git Status ==="
	git status --porcelain
	@echo
	@echo "=== Python Environment ==="
	$(PYTHON) --version
	@echo
	@echo "=== Installed Packages ==="
	$(PIP) list | head -20

# Debugging helpers
debug-build: ## Debug build issues
	$(PYTHON) -m build --wheel -v

debug-install: ## Debug installation issues
	$(PIP_INSTALL) -e . -v

# Performance testing
benchmark: ## Run performance benchmarks
	$(PYTHON) -m pytest $(EXAMPLES_DIR)/05_benchmark.ipynb --nbmake
