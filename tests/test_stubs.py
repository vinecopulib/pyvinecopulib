"""Tests for the `.pyi` stub generator (`scripts/generate_stubs.py`).

The generator is a build-time helper, not part of the installed package, so it
is loaded by file path. These tests guard the rendered output rather than the
gitignored stub artifacts on disk.
"""

import importlib.util
from pathlib import Path

from pyvinecopulib.families import BicopFamily

_SCRIPTS = (
  Path(__file__).resolve().parent.parent / "scripts" / "generate_stubs.py"
)


def _load_generator():
  spec = importlib.util.spec_from_file_location("generate_stubs", _SCRIPTS)
  assert spec is not None and spec.loader is not None
  mod = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(mod)
  return mod


def test_bicopfamily_members_rendered_as_class_attributes():
  """Enum members are declared inside the class body (see issue #223).

  ``pv.BicopFamily.clayton`` is the documented access pattern; the stub must
  declare each member as a typed class attribute so it passes static type
  checking, not only as a module-level constant.
  """
  gen = _load_generator()
  body = "\n".join(gen.render_class_stub(BicopFamily, "BicopFamily"))
  members = [
    n
    for n in dir(BicopFamily)
    if not n.startswith("_")
    and isinstance(getattr(BicopFamily, n), BicopFamily)
  ]
  assert members  # sanity: the enum exposes members
  for m in members:
    assert f"  {m}: BicopFamily = ..." in body
