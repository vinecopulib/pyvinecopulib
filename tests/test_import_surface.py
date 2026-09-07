"""Tests for what importing ``pyvinecopulib`` requires and exposes.

The package ships three optional subpackages behind extras. Which of them a
plain import, a star-import and `dir()` each reach is a contract users feel
directly -- a star-import that pulls in PyTorch fails on every machine that
does not have it -- and it is settled by one list, `__all__`, whose effect is
invisible from inside a fully-installed environment.

Blocking an extra with ``sys.modules[name] = None`` is how AGENTS.md says to
check this without a matrix run: the next ``import name`` raises, exactly as it
would where the package is absent.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

import pyvinecopulib as pv

#: The extras, and the subpackage each one gates.
_EXTRAS = ("torch", "sklearn", "scipy")


def _run(source: str) -> subprocess.CompletedProcess[str]:
  """Run ``source`` in a fresh interpreter, so no import is already cached."""
  return subprocess.run(
    [sys.executable, "-c", source], capture_output=True, text=True
  )


@pytest.mark.parametrize("extra", _EXTRAS)
def test_a_star_import_needs_no_optional_extra(extra: str) -> None:
  """`from pyvinecopulib import *` must not require an extra.

  It resolves every name in `__all__`, so listing a lazily-imported
  subpackage there quietly makes the whole extra a hard dependency of the one
  import form beginners reach for first.
  """
  done = _run(
    f"import sys; sys.modules[{extra!r}] = None\n"
    "exec('from pyvinecopulib import *')\n"
  )
  assert done.returncode == 0, done.stderr


@pytest.mark.parametrize("extra", _EXTRAS)
def test_a_plain_import_needs_no_optional_extra(extra: str) -> None:
  """The same for `import pyvinecopulib`, which has always held."""
  done = _run(
    f"import sys; sys.modules[{extra!r}] = None\nimport pyvinecopulib\n"
  )
  assert done.returncode == 0, done.stderr


def test_the_extras_are_absent_from_all_and_present_in_dir() -> None:
  """The two lists differ on purpose: `dir()` discovers, `__all__` binds."""
  for name in ("torch", "sklearn"):
    assert name not in pv.__all__
    assert name in dir(pv)


def test_the_extras_are_still_reachable_by_attribute() -> None:
  """Leaving them out of `__all__` must not make them unreachable."""
  pytest.importorskip("torch")
  assert pv.torch.__name__ == "pyvinecopulib.torch"
  pytest.importorskip("sklearn")
  assert pv.sklearn.__name__ == "pyvinecopulib.sklearn"


def test_margins_stays_in_all_because_it_needs_no_extra() -> None:
  """`margins` defers SciPy to the margin class that needs it."""
  assert "margins" in pv.__all__
  done = _run(
    "import sys; sys.modules['scipy'] = None\nimport pyvinecopulib.margins\n"
  )
  assert done.returncode == 0, done.stderr


def test_core_owns_the_margin_plumbing_and_margins_re_exports_it() -> None:
  """The layering: `core` owns the contract, `margins` the ecosystem adapters.

  `FitControlsMargin`, the coercion registry and the resolution helpers need no
  extra, and `core.VinedistBase` runs the two-step fit with them -- so keeping
  them a layer up had `core` importing *from* `margins`, privates included, at
  every use. They live in `core` and `margins` re-exports them, so the
  documented surface is unchanged.
  """
  import pyvinecopulib.core as core
  import pyvinecopulib.margins as margins

  for name in (
    "FitControlsMargin",
    "as_margin",
    "register_margin_adapter",
    "resolve_margins",
    "resolve_margin_controls",
  ):
    shared = getattr(margins, name)
    owner = getattr(shared, "__module__", "")
    assert owner.startswith("pyvinecopulib.core"), (name, owner)
  # And the one a caller reaches for beside its two siblings is right there.
  assert core.FitControlsMargin is margins.FitControlsMargin
  assert "FitControlsMargin" in dir(core)


def test_core_reaches_up_a_layer_only_where_it_must() -> None:
  """`core` must not import `margins` at module scope, and barely at all.

  Read statically rather than by watching `sys.modules`: importing
  `pyvinecopulib.core` runs the top-level `__init__`, which loads `margins`
  eagerly by design, so a runtime check would measure the wrong thing.

  Two deferred hops are expected and irreducible -- resolving the
  ``"parametric"`` alias and the ``"SciPyMargin"`` JSON ``kind`` to a class
  that lives behind an extra. Any more, or any at module scope, is the layer
  inversion coming back.
  """
  import ast
  import pathlib

  core = pathlib.Path("src/pyvinecopulib/core")
  if not core.is_dir():  # installed rather than checked out
    pytest.skip("source tree not available")

  module_scope: list[str] = []
  deferred: list[str] = []
  for path in sorted(core.glob("*.py")):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    top = set(tree.body)
    for node in ast.walk(tree):
      if not isinstance(node, (ast.Import, ast.ImportFrom)):
        continue
      target = getattr(node, "module", "") or ""
      if "margins" not in target:
        continue
      where = module_scope if node in top else deferred
      where.append(f"{path.name}:{node.lineno} -> {target}")

  assert module_scope == [], module_scope
  assert len(deferred) == 2, deferred
  # And neither reaches a private module of the layer above.
  assert not any("._" in d.split("-> ")[1] for d in deferred), deferred
