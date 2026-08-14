# Releasing pyvinecopulib

Maintainer checklist. This file is not part of the rendered user site.

Releases are tags on `main`. In the steady state the top heading of
`CHANGELOG.md` carries `(unreleased)`, so a released version is never
indistinguishable from an unreleased one.

## Checklist

1. **`main` is green** and every pull request intended for the release is
   merged.

2. **Open the release pull request.** In one commit:
   - retitle the top `CHANGELOG.md` heading from `(unreleased)` to the release
     date, e.g. `## 1.0.0 (2026-08-20)`;
   - set the version in `pyproject.toml` (`[project].version`), which is the
     single source of truth — CMake reads it through `SKBUILD_PROJECT_VERSION`
     and it reaches users as `pyvinecopulib.__version__`;
   - if `CITATION.cff` or a `version` field in `.zenodo.json` has been added
     by then, set the matching version there too.

3. **Check the metadata agrees**: `uv run python scripts/check_version.py`.
   CI runs the same check on every pull request.

4. **Merge**, then tag the merge commit:

   ```bash
   git checkout main && git pull
   git tag -a vX.Y.Z -m "pyvinecopulib X.Y.Z"
   git push origin vX.Y.Z
   ```

   > **The tag push publishes to PyPI, and PyPI never re-accepts a version
   > number.** A mistagged release cannot be undone — only superseded by a new
   > version. Confirm the tag points at the commit you mean before pushing.

5. **Confirm the release landed**: wheels and the sdist on PyPI, the GitHub
   release created from the changelog section, Read the Docs built the tag with
   `stable` repointed at it, and the Zenodo deposition.

6. **Open the next cycle immediately.** Add a fresh
   `## X.Y+1.0 (unreleased)` heading to `CHANGELOG.md` in a follow-up pull
   request, so the next change has somewhere to go.

7. **If the release surfaced a problem in the C++ library**, file it upstream
   rather than working around it here — the docstrings, signatures and CMake
   seams all belong to `vinecopulib`.

## What is automated

| Check | When |
|---|---|
| `scripts/check_version.py` | every pull request, via the `lint` job |
| wheels, sdist, tests, docs | every pull request |
| publish to PyPI + GitHub release | on a `v*` tag push |
| release drift (a dated changelog heading with no tag) | weekly |

## Version numbers

`pyvinecopulib` follows semantic versioning and tracks — but does not have to
match — the pinned `lib/vinecopulib` version. Honor the
[stability tiers](../AGENTS.md#stability-tiers): `core` / `families` / `utils`
prefer a deprecation alias over a hard break, while `sklearn` and `torch` may
break between minor releases.
