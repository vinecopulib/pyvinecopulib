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

4. **Rehearse the upload** (optional, but the only way to exercise the publish
   path before the irreversible one). Run the *Build Status* workflow manually
   from the Actions tab against the commit you are about to tag, with
   **Publish the built artifacts to Test PyPI** ticked. It builds the same
   wheels and sdist the tag would and uploads them to
   [Test PyPI](https://test.pypi.org/project/pyvinecopulib/).

   Test PyPI refuses a re-upload of a version it already has, exactly as PyPI
   does, so this is one shot per version — spend it on the commit you mean to
   tag, not on a work in progress. Needs the `test_pypi_password` secret.

5. **Run the mandatory CUDA rehearsal** on a GPU host before tagging. Run the
   CUDA device/compile tests and retain raw paired benchmark repetitions,
   confidence intervals, process RSS, and peak CUDA allocated/reserved memory
   with the release PR. The CPU matrix cannot validate these advertised paths.
   First verify `nvidia-smi`, then run
   `uv run pytest tests/test_torch_device.py tests/test_torch_bicop.py tests/test_torch_vinecop.py -m cuda --no-cov`
   and the maintained paired harness, `uv run python scripts/bench_torch_runner.py
   --repeats 10`; attach its raw JSON and the drivers' paired results to the PR.

6. **Merge**, then tag the merge commit:

   ```bash
   git checkout main && git pull
   git tag -a vX.Y.Z -m "pyvinecopulib X.Y.Z"
   git push origin vX.Y.Z
   ```

   > **The tag push publishes to PyPI, and PyPI never re-accepts a version
   > number.** A mistagged release cannot be undone — only superseded by a new
   > version. Confirm the tag points at the commit you mean before pushing.

7. **Confirm the release landed**: wheels and the sdist on PyPI, Read the Docs
   built the tag with `stable` repointed at it, and the Zenodo deposition. Then
   create the GitHub release manually from the matching changelog section; that
   release event is what Zenodo uses for its deposition when the repository is
   connected to Zenodo.

8. **Open the next cycle immediately.** Add a fresh
   `## X.Y+1.0 (unreleased)` heading to `CHANGELOG.md` in a follow-up pull
   request, so the next change has somewhere to go.

9. **If the release surfaced a problem in the C++ library**, file it upstream
   rather than working around it here — the docstrings, signatures and CMake
   seams all belong to `vinecopulib`.

## What is automated

| Check | When |
|---|---|
| `scripts/check_version.py` | every pull request, via the `lint` job |
| wheels, sdist, tests, docs | every pull request |
| `scripts/check_version.py --released --tag` | on a `v*` tag push, before publishing |
| publish to PyPI | on a `v*` tag push, after wheels, sdist, tests and docs |
| GitHub release + Zenodo deposition | maintainer action after PyPI publication |
| publish to Test PyPI | manually, from the Actions tab (see step 4) |

## Version numbers

`pyvinecopulib` follows semantic versioning and tracks — but does not have to
match — the pinned `lib/vinecopulib` version. Honor the
[stability tiers](../AGENTS.md#stability-tiers): `core` / `families` / `utils`
prefer a deprecation alias over a hard break, while `sklearn` and `torch` may
break between minor releases.
