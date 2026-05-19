# Thin wrappers around `uv`. See CONTRIBUTING.md.
.PHONY: help sync clean check format test test-examples docs sdist build notebooks
.DEFAULT_GOAL := help

UV := uv

help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-16s\033[0m %s\n", $$1, $$2}'

sync: ## Install all deps + editable build + pre-commit hooks
	$(UV) sync --all-extras --group dev --group test --group notebooks
	$(UV) pip install -e . --no-build-isolation
	$(UV) run pre-commit install

clean: ## Wipe build artifacts and Python caches
	rm -rf build/ dist/ *.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

check: ## Read-only lint + format-check + type-check (CI-safe)
	$(UV) run ruff check src tests
	$(UV) run ruff format --check src tests
	$(UV) run ty check

format: ## Apply ruff autofixes + format
	$(UV) run ruff check --fix src tests
	$(UV) run ruff format src tests

test: ## Run pytest suite
	$(UV) run pytest tests/

test-examples: ## Execute example notebooks as tests
	$(UV) run pytest --nbmake --nbmake-timeout=600 examples/

docs: ## Build HTML documentation
	$(UV) run sphinx-build -W -b html docs docs/_build/html

sdist: ## Build source distribution only
	$(UV) build --sdist

build: ## Build sdist and wheel
	$(UV) build

notebooks: ## Re-execute example notebooks
	$(UV) run python scripts/regenerate_notebooks.py
