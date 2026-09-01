# Thin wrappers around `uv`. See CONTRIBUTING.md.
.PHONY: help sync clean lint check format test test-examples docs sdist build notebooks
.DEFAULT_GOAL := help

UV := uv
# Every `uv run` otherwise re-syncs, which reinstalls the project *with*
# build isolation and so re-poisons `build/` (see `sync`). The workflow is
# `make sync` first; CI already sets this for the targets it runs.
export UV_NO_SYNC := 1

help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-16s\033[0m %s\n", $$1, $$2}'

sync: ## Install all deps + editable build + pre-commit hooks
# `--no-install-project`: `build-dir` is a persistent tree, so a first,
# build-isolated install of the project bakes that throwaway environment's
# `ninja` path into `CMakeCache.txt`. The editable rebuild then invokes a
# binary that no longer exists, and every later import dies in cmake.
	$(UV) sync --no-install-project --all-extras --group dev --group test --group notebooks
	$(UV) pip install -e . --no-build-isolation --python .venv
	$(UV) run pre-commit install

# The staged docs sources are wiped too. `conf.py` stages repo files into the
# docs tree and writes two toctrees beside them, and it skips anything already
# there -- so a copy left by another branch is never refreshed, and the
# nitpicky build fails on a file this branch does not have.
clean: ## Wipe build artifacts, staged docs sources, and Python caches
	rm -rf build/ dist/ *.egg-info/
	rm -rf docs/_build docs/_generate docs/examples docs/examples.rst \
	       docs/features.rst docs/README.md docs/_README_inlined.md \
	       docs/CHANGELOG.md
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

lint: ## Lint + format-check + security-lint; needs no compiled extension
	$(UV) run ruff check src tests
	$(UV) run ruff format --check src tests
	$(UV) run bandit -c pyproject.toml -q -r src/pyvinecopulib scripts
	$(UV) run codespell src tests scripts docs examples .github *.md *.cff

check: lint ## `lint` plus the type check, which reads the generated .pyi stubs
	$(UV) run ty check

format: ## Apply ruff autofixes + format
	$(UV) run ruff check --fix src tests
	$(UV) run ruff format src tests

test: ## Run pytest suite
	$(UV) run pytest tests/

test-examples: ## Execute example notebooks as tests
	@command -v dot >/dev/null || { \
	  echo "Graphviz 'dot' is required for notebook execution; install Graphviz first."; \
	  exit 2; \
	}
	PYTHONWARNINGS="error::DeprecationWarning:__main__" \
	  $(UV) run pytest --nbmake --nbmake-timeout=600 examples/

docs: ## Build HTML documentation
	$(UV) run sphinx-build -W -b html docs docs/_build/html

sdist: ## Build source distribution only
	$(UV) build --sdist

build: ## Build sdist and wheel
	$(UV) build

notebooks: ## Re-execute example notebooks
	$(UV) run python scripts/regenerate_notebooks.py
