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
