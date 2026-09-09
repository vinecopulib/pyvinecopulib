"""Tests for the `.pyi` stub generator (`scripts/generate_stubs.py`).

The generator is a build-time helper, not part of the installed package, so it
is loaded by file path. These tests guard the rendered output rather than the
gitignored stub artifacts on disk.
"""

import builtins
import ast
import importlib.util
from pathlib import Path
from types import ModuleType
from typing import cast

import pyvinecopulib as pv
from pyvinecopulib.families import BicopFamily

_SCRIPTS = (
  Path(__file__).resolve().parent.parent / "scripts" / "generate_stubs.py"
)


def _load_generator() -> ModuleType:
  spec = importlib.util.spec_from_file_location("generate_stubs", _SCRIPTS)
  assert spec is not None and spec.loader is not None
  mod = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(mod)
  return mod


def test_bicopfamily_members_rendered_as_class_attributes() -> None:
  """Enum members are declared inside the class body (see issue #223).

  ``pv.BicopFamily.clayton`` is the documented access pattern; the stub must
  declare each member as a typed class attribute so it passes static type
  checking, not only as a module-level constant.
  """
  gen = _load_generator()
  body = "\n".join(
    cast("list[str]", gen.render_class_stub(BicopFamily, "BicopFamily"))
  )
  members = [
    n
    for n in dir(BicopFamily)
    if not n.startswith("_")
    and isinstance(getattr(BicopFamily, n), BicopFamily)
  ]
  assert members  # sanity: the enum exposes members
  for m in members:
    assert f"  {m}: BicopFamily = ..." in body


def _rendered(cls: type, name: str) -> str:
  gen = _load_generator()
  # The generator is loaded by file path, so its return type is opaque here.
  return "\n".join(cast("list[str]", gen.render_class_stub(cls, name)))


def test_keyword_only_arguments_survive_into_the_stub() -> None:
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
    (bicop, "def sample(self, n: int | None = None"),
    (bicop, "*, parameters:"),
    (vinecop, "def rosenblatt(self"),
    (vinecop, "def inverse_rosenblatt(self"),
    (vinecop, "def sample_conditional(self"),
  ]:
    assert needle in haystack, f"missing {needle!r}"

  # `conditioning_set` must appear after a bare `*` on all three methods, so
  # position 2 keeps meaning `num_threads` / `qrng`.
  for method in ("rosenblatt", "inverse_rosenblatt", "sample_conditional"):
    line = next(
      line for line in vinecop.splitlines() if f"def {method}(self" in line
    )
    assert "*, conditioning_set:" in line, line


def test_a_deprecated_alias_renders_with_a_real_signature() -> None:
  """A shim must not degrade to ``*args, **kwargs`` in the stubs.

  `generate_stubs.py` recovers a signature -- and ``@staticmethod`` -- by parsing
  the first line of the docstring, so the alias copies that line from the method
  it forwards to. If that ever stops working, a type checker sees an untyped
  callable and `RVineStructure.simulate` loses its staticmethod, which is exactly
  what this asserts against.
  """
  from pyvinecopulib.core import Bicop, RVineStructure

  bicop = _rendered(Bicop, "Bicop")
  structure = _rendered(RVineStructure, "RVineStructure")

  assert "def simulate(self, n: int | None = None" in bicop
  assert "*args" not in bicop.split("def simulate")[1].split("\n")[0]

  # The static alias keeps its decorator and its `self`-less signature.
  line = next(
    line for line in structure.splitlines() if "def simulate(" in line
  )
  assert line.lstrip().startswith("def simulate(d: int"), line


def test_from_data_accepts_a_dynamically_sized_matrix() -> None:
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


def test_no_binding_is_an_overload_set() -> None:
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


def test_subclass_declares_its_base() -> None:
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


def test_base_without_a_definition_here_is_not_declared() -> None:
  """Only bases the stub itself defines are named.

  ``BicopFamily`` derives from ``enum.Enum``, which the stub never declares, so
  naming it would leave a dangling reference.
  """
  gen = _load_generator()
  body = "\n".join(
    cast("list[str]", gen.render_class_stub(BicopFamily, "BicopFamily"))
  )
  assert body.startswith("class BicopFamily:")


def test_docstrings_are_escaped_as_valid_python() -> None:
  """LaTeX backslashes cannot form truncated escapes in generated stubs."""
  gen = _load_generator()
  rendered = "\n".join(gen.render_docstring(r"criterion $\xi$", 0))
  ast.parse(rendered)
  assert r"\\xi" in rendered


def test_every_generated_stub_parses() -> None:
  """The build's PEP 561 artifacts are syntactically valid Python."""
  root = Path(pv.__file__).resolve().parent
  stubs = sorted(root.glob("**/__init__.pyi"))
  assert stubs
  for stub in stubs:
    ast.parse(stub.read_text(encoding="utf-8"), filename=str(stub))


def test_every_stub_annotation_resolves() -> None:
  """A name an annotation refers to has to be defined or imported.

  Parsing does not catch this: `~ArrayT` is a valid unary invert on a name and
  `Optional[module]` is a valid subscript, so a stub that renders
  `repr(TypeVar)` or a class's bare `__name__` parses cleanly and then means
  nothing to a consumer. Both shapes reached the published `core` stub, and
  the second silently degraded every generic signature in it to `Any`.
  """
  root = Path(pv.__file__).resolve().parent
  stubs = sorted(root.glob("**/__init__.pyi"))
  assert stubs

  unresolved: dict[str, set[str]] = {}
  contributed: dict[str, int] = {}
  for stub in stubs:
    text = stub.read_text(encoding="utf-8")
    tree = ast.parse(text, filename=str(stub))

    bound: set[str] = set(dir(builtins)) | {"typing"}
    for node in ast.walk(tree):
      if isinstance(node, ast.Import):
        bound |= {(a.asname or a.name).split(".")[0] for a in node.names}
      elif isinstance(node, ast.ImportFrom):
        bound |= {a.asname or a.name for a in node.names}
      elif isinstance(node, (ast.ClassDef, ast.FunctionDef)):
        bound.add(node.name)
      elif isinstance(node, ast.AnnAssign) and isinstance(
        node.target, ast.Name
      ):
        bound.add(node.target.id)
      elif isinstance(node, ast.Assign):
        bound |= {t.id for t in node.targets if isinstance(t, ast.Name)}

    # Every annotation, whether written as a string or evaluated into one.
    referenced: set[str] = set()
    for node in ast.walk(tree):
      annotations = []
      if isinstance(node, ast.arg) and node.annotation is not None:
        annotations.append(node.annotation)
      elif isinstance(node, ast.FunctionDef) and node.returns is not None:
        annotations.append(node.returns)
      elif isinstance(node, ast.AnnAssign):
        annotations.append(node.annotation)
      for annotation in annotations:
        source = (
          annotation.value
          if isinstance(annotation, ast.Constant)
          and isinstance(annotation.value, str)
          else ast.unparse(annotation)
        )
        try:
          inner = ast.parse(source, mode="eval")
        except SyntaxError:
          unresolved.setdefault(stub.name, set()).add(source)
          continue
        referenced |= {n.id for n in ast.walk(inner) if isinstance(n, ast.Name)}

    missing = {n for n in referenced if n not in bound}
    if missing:
      unresolved.setdefault(stub.name, set()).update(missing)
    # An empty `missing` is also what a walk that collected nothing reports,
    # so every stub carrying a signature has to have contributed a reference
    # before the empty result means anything.
    if "def " in text:
      contributed[str(stub.relative_to(root))] = len(referenced)

  assert unresolved == {}, {k: sorted(v) for k, v in unresolved.items()}
  assert contributed, "no stub carried a signature: the glob found nothing"
  barren = sorted(k for k, n in contributed.items() if n == 0)
  assert barren == [], f"annotations collected from no names in: {barren}"
