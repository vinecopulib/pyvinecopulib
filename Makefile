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
# build-isolated install of the project writes that throwaway environment's
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

#: Everything whose Python this repository owns. Not `.`: see `lint`.
OWNED := src tests examples scripts docs

lint: ## Lint + format-check + security-lint; needs no compiled extension
# Both linters name what this repository owns rather than scanning the tree:
# CI's dependency step unpacks Eigen and Boost into the workspace, and Eigen
# ships Python 2 scripts that ruff reads as syntax errors.
	$(UV) run ruff check $(OWNED)
	$(UV) run ruff format --check $(OWNED)
	$(UV) run bandit -c pyproject.toml -q -r src/pyvinecopulib scripts
# Tracked files, for the same reason plus one more: a shell glob would also
# read whatever untracked notes are in the working tree, so a local run could
# fail where CI -- which checks out only tracked files -- passes.
	git ls-files -z ':!:lib/*' | xargs -0 $(UV) run codespell --

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
# autosummary adds a stub per documented name under `_generate` and never
# removes one, so a name that stops being public leaves an orphan stub whose
# `autofunction` then fails the `-W` build. Both paths are generated, and
# `_build` is left alone so the HTML stays incremental.
	rm -rf docs/_generate docs/features.rst
	$(UV) run sphinx-build -W -b html docs docs/_build/html

sdist: ## Build source distribution only
	$(UV) build --sdist

build: ## Build sdist and wheel
	$(UV) build

notebooks: ## Re-execute example notebooks
	$(UV) run python scripts/regenerate_notebooks.py
