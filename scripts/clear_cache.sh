#!/bin/zsh

setopt +o nomatch  # allow unmatched globs

RM=/bin/rm  # fallback to full path

# List of paths to remove
paths=(
  build/
  _skbuild/
  dist/
  .pytest_cache/
  .mypy_cache/
  .ruff_cache/
  __pycache__/
  src/**/__pycache__/
  scripts/**/__pycache__/
  docs/**/__pycache__/
  tests/**/__pycache__/
  pyvinecopulib/__pycache__/
  CMakeCache.txt
  CMakeFiles/
  .cache/
  *.egg-info
)

for path in $paths; do
  $RM -rf $path
done

echo "✅ Cleanup complete."