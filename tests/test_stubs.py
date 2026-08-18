"""Tests for the `.pyi` stub generator (`scripts/generate_stubs.py`).

The generator is a build-time helper, not part of the installed package, so it
is loaded by file path. These tests guard the rendered output rather than the
gitignored stub artifacts on disk.
"""

import importlib.util
from pathlib import Path

import pyvinecopulib as pv
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


def _rendered(cls, name: str) -> str:
  gen = _load_generator()
  return "\n".join(gen.render_class_stub(cls, name))


def test_keyword_only_arguments_survive_into_the_stub():
  """The new conditioning / per-row arguments are keyword-only.

  A nanobind overload set renders as ``"Overloaded function."`` and loses its
  signature entirely, so each of these is bound as one method with an internal
  dispatch. Pin the rendered signature: a dropped ``*`` or a renamed
  ``nb::arg`` would otherwise reach users' type checkers unnoticed.
  """
  from pyvinecopulib.core import Bicop, Vinecop

  bicop = _rendered(Bicop, "Bicop")
  vinecop = _rendered(Vinecop, "Vinecop")

  for haystack, needle in [
    (bicop, "def simulate(self, n: int | None = None"),
    (bicop, "*, parameters:"),
    (vinecop, "def rosenblatt(self"),
    (vinecop, "def inverse_rosenblatt(self"),
    (vinecop, "def simulate_conditional(self"),
  ]:
    assert needle in haystack, f"missing {needle!r}"

  # `conditioning_set` must appear after a bare `*` on all three methods, so
  # position 2 keeps meaning `num_threads` / `qrng`.
  for method in ("rosenblatt", "inverse_rosenblatt", "simulate_conditional"):
    line = next(
      line for line in vinecop.splitlines() if f"def {method}(self" in line
    )
    assert "*, conditioning_set:" in line, line


def test_from_data_accepts_a_dynamically_sized_matrix():
  """`Bicop.from_data` must not re-acquire a static two-column shape.

  A statically two-column Eigen type makes the discrete layouts unpassable.
  """
  from pyvinecopulib.core import Bicop

  line = next(
    line
    for line in _rendered(Bicop, "Bicop").splitlines()
    if "def from_data(" in line
  )
  assert "shape=(*, 2)" not in line, line


def test_no_binding_is_an_overload_set():
  """Alternative constructors are named factories, not C++-style overloads.

  An overload set costs twice: nanobind concatenates the docstrings, so two
  numpydoc ``Parameters`` sections collide and fail the docs build, and the
  generator above renders only the first signature, so a type checker rejects
  every call matching the others. See the convention in AGENTS.md.
  """
  import inspect
  import re

  from pyvinecopulib import core, families, utils

  overloaded = []
  for module in (core, families, utils):
    for name in getattr(module, "__all__", []):
      obj = getattr(module, name)
      if not inspect.isclass(obj):
        continue
      for attr_name in dir(obj):
        try:
          attr = getattr(obj, attr_name)
        except Exception:
          continue
        doc = inspect.getdoc(attr) or ""
        pattern = re.compile(rf"^{re.escape(attr_name)}\(.*\)\s*->")
        signatures = 0
        for line in doc.splitlines():
          if not pattern.match(line.strip()):
            break
          signatures += 1
        if signatures > 1:
          overloaded.append(f"{module.__name__}.{name}.{attr_name}")

  assert not overloaded, (
    "bound as overload sets; bind a named factory instead: "
    + ", ".join(sorted(set(overloaded)))
  )


def test_subclass_declares_its_base():
  """``DVineStructure`` / ``CVineStructure`` derive from ``RVineStructure``.

  The binding declares the inheritance and the runtime MRO carries it, so the
  stub has to as well: passing a D-vine wherever an ``RVineStructure`` is
  expected is the documented way to use one, and a stub that omits the base
  makes that a type error.
  """
  gen = _load_generator()
  for cls in (pv.DVineStructure, pv.CVineStructure):
    body = "\n".join(gen.render_class_stub(cls, cls.__name__))
    assert body.startswith(f"class {cls.__name__}(RVineStructure):")


def test_base_without_a_definition_here_is_not_declared():
  """Only bases the stub itself defines are named.

  ``BicopFamily`` derives from ``enum.Enum``, which the stub never declares, so
  naming it would leave a dangling reference.
  """
  gen = _load_generator()
  body = "\n".join(gen.render_class_stub(BicopFamily, "BicopFamily"))
  assert body.startswith("class BicopFamily:")
