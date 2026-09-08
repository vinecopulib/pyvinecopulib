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

import ast
import pathlib
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

  Three deferred hops are expected and irreducible, each resolving a name to a
  class that lives behind an extra: the ``"parametric"`` alias, the
  ``"SciPyMargin"`` JSON ``kind``, and the adapter ``as_margin`` builds for an
  OpenTURNS object. Any more, or any at module scope, is the layer inversion
  coming back.
  """
  import ast
  import pathlib

  core = pathlib.Path("src/pyvinecopulib/core")
  if not core.is_dir():  # installed rather than checked out
    pytest.skip("source tree not available")

  def targets_of(node: ast.Import | ast.ImportFrom) -> list[str]:
    """Resolve one import statement to the dotted modules it names.

    Resolved rather than matched as a substring, because ``core`` has a
    sibling named ``_margins``: level 1 is a module inside ``core``, level 2
    the package above it.
    """
    if isinstance(node, ast.Import):
      return [alias.name for alias in node.names]
    here = ["pyvinecopulib", "core"][: 2 - ((node.level or 1) - 1)]
    return [".".join([*here, node.module or ""])]

  module_scope: list[str] = []
  deferred: list[str] = []
  for path in sorted(core.glob("*.py")):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    top = set(tree.body)
    for node in ast.walk(tree):
      if not isinstance(node, (ast.Import, ast.ImportFrom)):
        continue
      for target in targets_of(node):
        if not target.startswith("pyvinecopulib.margins"):
          continue
        where = module_scope if node in top else deferred
        where.append(f"{path.name}:{node.lineno} -> {target}")

  assert module_scope == [], module_scope
  assert len(deferred) == 3, deferred
  # And neither reaches a private module of the layer above.
  assert not any("._" in d.split("-> ")[1] for d in deferred), deferred


def test_the_fit_callback_aliases_resolve_from_any_module() -> None:
  """`FitEdge` / `FitLevel` must carry no unresolved forward reference.

  A string inside a type *alias* is resolved in whichever module uses the
  alias, not where it is defined -- so quoting one made every importer keep a
  name in scope for it, and the docs build, which resolves annotations at
  runtime, was what noticed when one stopped.
  """
  import typing

  from pyvinecopulib.core.vinecop_base import VinecopBase

  for method in (VinecopBase.fit, VinecopBase.select, VinecopBase.from_data):
    typing.get_type_hints(method)


def test_agents_md_names_every_module_that_exists() -> None:
  """The package tree in AGENTS.md is a map, so it has to match the territory.

  It was missing six modules -- including `_validation.py` and `_trim.py`, the
  layout and domain steps of the input pipeline the file spends a page on.
  Comparing the two by eye is what let that happen.
  """
  import os
  import pathlib
  import re

  spec = pathlib.Path("AGENTS.md")
  root = pathlib.Path("src/pyvinecopulib")
  if not spec.is_file() or not root.is_dir():
    pytest.skip("source tree not available")

  tree = re.search(r"```text\n(.*?)```", spec.read_text(encoding="utf-8"), re.S)
  assert tree is not None, "AGENTS.md has no package-structure block"
  listed = {
    os.path.basename(token)
    for line in tree.group(1).splitlines()
    for token in re.findall(r"[\w./*]+\.py", line.split("#")[0])
  }
  actual = {
    path.name
    for path in root.rglob("*.py")
    if path.name != "__init__.py" and "__pycache__" not in path.parts
  }
  assert actual - listed == set(), sorted(actual - listed)


def test_agents_md_public_api_lists_match_the_code() -> None:
  """AGENTS.md's "Public APIs" section is a second copy of every `__all__`.

  Two copies drift, and this one had: it placed `Kde1d` in `utils`, where it
  has never been.
  """
  import pathlib
  import re

  spec = pathlib.Path("AGENTS.md")
  if not spec.is_file():
    pytest.skip("source tree not available")
  text = spec.read_text(encoding="utf-8")

  import pyvinecopulib.core as core
  import pyvinecopulib.families as families
  import pyvinecopulib.margins as margins
  import pyvinecopulib.utils as utils

  subpackages = {"core", "families", "utils", "margins", "sklearn", "torch"}
  section = text[text.index("## Public APIs") :]
  for label, module in (
    ("`pyvinecopulib.core`", core),
    ("`pyvinecopulib.families`", families),
    ("`pyvinecopulib.utils`", utils),
    ("`pyvinecopulib.margins`", margins),
  ):
    start = section.index(f"- **{label}**")
    entry = section[start : section.index("\n- **", start + 1)]
    named = set(re.findall(r"`([A-Za-z_][A-Za-z_0-9]*)`", entry))
    exported = {
      name
      for name in getattr(module, "__all__", ())
      if not name.startswith("__")
    } - subpackages
    # Both directions: a name the code exports and the list omits is the drift
    # that hid `Kde1d`, and a name the list claims and the code does not export
    # is the drift that put `Kde1d` in `utils`. One-way would catch neither.
    assert exported - named == set(), (
      label,
      "undocumented",
      sorted(exported - named),
    )
    stale = {n for n in named if n in dir(module)} | (named & exported)
    invented = {
      n
      for n in named
      if n not in exported
      and n not in dir(module)
      and n[:1].isupper()
      and n not in subpackages
    }
    assert invented == set(), (label, "claimed but absent", sorted(invented))
    del stale


#: The package root, which sits above every layer.
_TOP = "<top>"

#: Every intra-package import edge, and whether every occurrence of it is
#: function-local. An edge missing here is a layer crossing nobody argued
#: for; an edge whose value flips from `True` is a deferral that stopped
#: deferring, which is how an optional extra becomes a hard dependency.
_LAYER_EDGES: dict[tuple[str, str], bool] = {
  # The root re-export surface, above everything.
  (_TOP, "core"): False,
  (_TOP, "families"): False,
  (_TOP, "margins"): False,
  (_TOP, "utils"): False,
  (_TOP, "pyvinecopulib_ext"): False,
  (_TOP, "_cpu"): False,
  (_TOP, "_deprecations"): False,
  # Tier 2 -> tier 1, plus the two edges within tier 2. Constructing
  # `TorchVinecopBackend` is the opt-in signal that PyTorch is required, so
  # that one import has to stay inside the constructor.
  ("margins", "core"): False,
  ("torch", "core"): False,
  ("torch", "utils"): False,
  ("torch", "pyvinecopulib_ext"): False,
  ("sklearn", "core"): False,
  ("sklearn", "margins"): False,
  ("sklearn", "torch"): True,
  # Tier 1, and its two deferred hops: up into `margins` for `SciPyMargin`
  # (see `test_core_reaches_up_a_layer_only_where_it_must`), and across into
  # the binding's helpers.
  ("core", "pyvinecopulib_ext"): False,
  ("core", "_deprecations"): False,
  ("core", "_python_helpers"): True,
  ("core", "margins"): True,
  ("families", "pyvinecopulib_ext"): False,
  ("utils", "pyvinecopulib_ext"): False,
  ("utils", "_python_helpers"): False,
  ("_python_helpers", "core"): False,
  # The x86-64-v3 guard reads the build's own record of what it compiled.
  ("_cpu", "_build_info"): True,
}


def _measure_layer_edges(root: pathlib.Path) -> dict[tuple[str, str], bool]:
  """Read every intra-package import edge, and whether it is deferred.

  Parameters
  ----------
  root : pathlib.Path
      The ``src/pyvinecopulib`` directory.

  Returns
  -------
  dict
      ``(importer, imported) -> whether every occurrence is function-local``.
  """
  layers = {p.name for p in root.iterdir() if (p / "__init__.py").is_file()}
  layers |= {p.stem for p in root.glob("*.py") if p.stem != "__init__"}
  # Two layers the build produces rather than the checkout: the extension
  # module is never a `.py`, and `_build_info.py` is generated. Both are
  # layers whether or not this tree has been built, or the edges into them
  # would read as absent on a fresh clone.
  layers |= {"pyvinecopulib_ext", "_build_info"}

  def _layer(parts: tuple[str, ...]) -> str:
    if parts == ("__init__.py",):
      return _TOP
    return parts[0] if parts[0] in layers else parts[0].removesuffix(".py")

  def _targets(
    node: ast.Import | ast.ImportFrom, parts: tuple[str, ...]
  ) -> set[str]:
    """The layers one import statement names, whatever form it takes."""
    names = [a.name for a in node.names]
    if isinstance(node, ast.Import):
      # `import pyvinecopulib.torch` -- every alias, not just the first.
      return {n.split(".")[1] for n in names if n.startswith("pyvinecopulib.")}
    if node.level:  # `from .. import x`, `from ..core import y`
      base = list(parts[:-1])
      if node.level > 1:
        base = base[: len(base) - (node.level - 1)]
      prefix = base + (node.module.split(".") if node.module else [])
    elif (absolute := node.module or "").startswith("pyvinecopulib"):
      prefix = absolute.split(".")[1:]
    else:
      return set()
    # With a prefix, the first component is the layer; without one the
    # statement is `from . import <layer>`, so the names are the layers.
    return {prefix[0]} if prefix else set(names)

  edges: dict[tuple[str, str], bool] = {}
  for path in sorted(root.rglob("*.py")):
    parts = path.relative_to(root).parts
    src = _layer(parts)
    if src not in layers | {_TOP}:
      continue
    tree = ast.parse(path.read_text(encoding="utf-8"))
    top = set(tree.body)
    for node in ast.walk(tree):
      if not isinstance(node, (ast.Import, ast.ImportFrom)):
        continue
      for dst in _targets(node, parts) & layers:
        if dst == src:
          continue
        key = (src, dst)
        edges[key] = edges.get(key, True) and node not in top
  return edges


def test_the_layers_only_depend_downwards() -> None:
  """The dependency direction AGENTS.md draws is the one the code has.

  A new import is the cheapest way to invert a layer, and the damage shows up
  far away -- as an extra a plain `import pyvinecopulib` suddenly needs. So
  the whole edge set is pinned, not just the one edge that inverted once:
  adding an import that crosses layers has to be an edit here, with
  the reason written beside it.
  """
  root = pathlib.Path("src/pyvinecopulib")
  if not root.is_dir():  # installed rather than checked out
    pytest.skip("source tree not available")

  measured = _measure_layer_edges(root)
  assert set(measured) == set(_LAYER_EDGES), {
    "undeclared": sorted(set(measured) - set(_LAYER_EDGES)),
    "gone": sorted(set(_LAYER_EDGES) - set(measured)),
  }
  must_defer = {e for e, deferred in _LAYER_EDGES.items() if deferred}
  eager = {e for e in must_defer if not measured[e]}
  assert eager == set(), eager


def test_no_test_guards_an_extra_by_a_first_party_import() -> None:
  """`importorskip` must name the extra, never a `pyvinecopulib` module.

  Whether importing one of our own modules raises depends on its internals --
  `pyvinecopulib.margins` imports fine without SciPy, since only
  `SciPyMargin`'s constructor needs it -- so such a guard skips nothing and
  the test then fails on the extras-free CI legs instead. Naming the extra
  cannot go stale that way.
  """
  offenders: list[str] = []
  for path in sorted(pathlib.Path("tests").glob("*.py")):
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
      if not isinstance(node, ast.Call):
        continue
      func = node.func
      named = getattr(func, "attr", None) or getattr(func, "id", None)
      if named != "importorskip" or not node.args:
        continue
      first = node.args[0]
      if (
        isinstance(first, ast.Constant)
        and isinstance(first.value, str)
        and first.value.startswith("pyvinecopulib")
      ):
        offenders.append(f"{path.name}:{node.lineno} -> {first.value}")
  assert offenders == [], offenders
