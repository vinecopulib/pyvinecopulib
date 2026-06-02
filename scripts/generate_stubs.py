#!/usr/bin/env python3
"""Generate `.pyi` stubs for pyvinecopulib subpackages.

CMake POST_BUILD invokes this after linking the extension. Sources and the
freshly built `.so` are staged into a tempdir so importlib sees a complete
package. Modules whose import raises ImportError (e.g. `pyvinecopulib.sklearn`
without scikit-learn) are skipped with a warning.
"""

import argparse
import importlib
import inspect
import re
import shutil
import sys
import tempfile
from pathlib import Path
from types import BuiltinFunctionType, FunctionType
from typing import Optional


def wrap_known_types(sig: str, known_types: set[str]) -> str:
  for name in known_types:
    sig = re.sub(rf"= *\[{name}\.\w+(?: *, *{name}\.\w+)*\]", "= ...", sig)
    sig = re.sub(rf"= *{name}\s*\(\)", "= ...", sig)
    sig = re.sub(rf"= *{name}\s*\.\w+", "= ...", sig)
    sig = re.sub(rf"\b{name}\b(?!\s*\()", f'"{name}"', sig)
  return sig


def render_python_function_stub(
  fct, name: str, known_types: Optional[set[str]] = None, indent: int = 2
) -> list[str]:
  """Render stub for regular Python functions using inspect.signature()"""
  import inspect

  lines = []
  indent_str = " " * indent

  try:
    sig = inspect.signature(fct)
    sig_str = str(sig)
    if known_types:
      sig_str = wrap_known_types(sig_str, known_types)
    lines.append(f"def {name}{sig_str}:")
  except Exception:
    lines.append(f"def {name}(*args, **kwargs):")

  doc = inspect.getdoc(fct)
  if doc:
    lines.append(f'{indent_str}"""')
    for line in doc.splitlines():
      lines.append(f"{indent_str}{line}")
    lines.append(f'{indent_str}"""')

  lines.append(f"{indent_str}...\n")
  return lines


def render_nanobind_function_stub(
  fct, name: str, known_types: Optional[set[str]] = None, indent: int = 2
) -> list[str]:
  import inspect
  import re

  doc = inspect.getdoc(fct) or ""
  lines = []
  doc_lines = doc.splitlines()
  first_line = doc_lines[0].strip() if doc_lines else ""

  sig_pattern = re.compile(rf"^{re.escape(name)}\((.*?)\)\s*->\s*(.*)$")
  match = sig_pattern.match(first_line)

  if match:
    args, ret = match.groups()
    if known_types:
      args = wrap_known_types(args, known_types)
      ret = wrap_known_types(ret, known_types)
    lines.append(f"def {name}({args}) -> {ret}:")
    remaining_doc = doc_lines[1:]
  else:
    lines.append(f"def {name}(*args, **kwargs):")
    remaining_doc = doc_lines

  indent_str = " " * indent
  if remaining_doc:
    lines.append(f'{indent_str}"""')
    for line in remaining_doc:
      lines.append(f"{indent_str}{line}")
    lines.append(f'{indent_str}"""')

  lines.append(f"{indent_str}...\n")
  return lines


def infer_method_decorator(name: str, docstring: str) -> Optional[str]:
  """Infer whether a method should be decorated as static or classmethod."""
  if not docstring:
    return []

  first_line = docstring.strip().splitlines()[0]
  sig_pattern = re.compile(rf"^{re.escape(name)}\((.*?)\)\s*(->.*)?$")
  match = sig_pattern.match(first_line)
  if not match:
    return []

  args_str = match.group(1)
  if not args_str.strip():
    return []

  args = []
  for arg in args_str.split(","):
    name_part = arg.strip().split(":", 1)[0].strip()
    if name_part:
      args.append(name_part)
  if not args:
    return "@staticmethod"
  if "cls" in args and "self" not in args:
    return "@classmethod"
  if "self" not in args:
    return "@staticmethod"
  return None  # regular instance method


def render_class_stub(
  cls, name: str, known_types: Optional[set[str]] = None, indent: int = 2
) -> list[str]:
  lines = [f"class {name}:"]
  doc = inspect.getdoc(cls)
  inner_indent = " " * indent
  if doc:
    lines.append(f'{inner_indent}"""')
    for line in doc.splitlines():
      lines.append(f"{inner_indent}{line}")
    lines.append(f'{inner_indent}"""')
  else:
    lines.append(f"{inner_indent}...")
    return lines

  for attr_name in dir(cls):
    if attr_name.startswith("_") and attr_name != "__init__":
      continue

    try:
      attr = getattr(cls, attr_name)
    except Exception as exc:
      # Some bound attributes raise on access (e.g. nanobind descriptors);
      # skip them but log so the omission isn't silent.
      print(f"skip {cls.__name__}.{attr_name}: {exc}", file=sys.stderr)
      continue

    if callable(attr) and hasattr(attr, "__doc__"):
      doc = inspect.getdoc(attr) or ""
      decorator = infer_method_decorator(attr_name, doc)
      if decorator:
        lines.append(f"\n{inner_indent}{decorator}")
      fct_lines = [
        inner_indent + line
        for line in render_nanobind_function_stub(
          attr, attr_name, known_types=known_types, indent=indent
        )
      ]
      lines.extend(fct_lines)
    elif isinstance(attr, property) or inspect.isdatadescriptor(attr):
      lines.append(f"{inner_indent}@property")
      doc = inspect.getdoc(attr) or ""
      lines.append(f"{inner_indent}def {attr_name}(self) -> Any: ...")
      if attr.fset is not None:
        lines.append(f"{inner_indent}@{attr_name}.setter")
        lines.append(
          f"{inner_indent}def {attr_name}(self, value: Any) -> None: ..."
        )

  lines.append("")  # Final newline
  return lines


def cleanup_stub(stub: str) -> str:
  # In `def …(...) -> ...:` lines: inputs become `ArrayLike` (lists, scalars,
  # buffers), returns become `NDArray[Any]` (subscriptable, supports arithmetic).
  # Match `numpy.ndarray[ … ]` / `np.ndarray[ … ]` with up to two levels of
  # nested brackets in the parameter list (covers `ndarray[tuple[Any, ...],
  # dtype[Any]]`-style annotations from `inspect.signature` / `NDArray[Any]`).
  ndarray_pattern = re.compile(
    r"(?:numpy|np)\.ndarray\[(?:[^\[\]]|\[(?:[^\[\]]|\[[^\[\]]*\])*\])*\]"
    r"|\bnumpy\.ndarray\b|\bnp\.ndarray\b"
  )

  def _rewrite_def_line(match: re.Match) -> str:
    signature, return_annot = match.group(1), match.group(2)
    signature = ndarray_pattern.sub("ArrayLike", signature)
    return_annot = ndarray_pattern.sub("NDArray[Any]", return_annot)
    return f"def {signature} -> {return_annot}:"

  stub = re.sub(r"def ([^\n]+?) -> ([^\n]+?):", _rewrite_def_line, stub)
  stub = ndarray_pattern.sub("ArrayLike", stub)

  stub = re.sub(
    r"Union\[numpy\._typing\._array_like\._Buffer, numpy\._typing\._array_like\._SupportsArray\[numpy\.dtype\[Any\]\], numpy\._typing\._nested_sequence\._NestedSequence\[numpy\._typing\._array_like\._SupportsArray\[numpy\.dtype\[Any\]\]\], complex, bytes, str, numpy\._typing\._nested_sequence\._NestedSequence\[complex \| bytes \| str\]\]",
    "ArrayLike",
    stub,
  )

  stub = re.sub(r"= *(?:np\.)?array\(\[\], *dtype=.*?\)", "= ...", stub)
  stub = re.sub(
    r"= *(?:np\.)?array\(\[\], *(?:shape=\([^\)]*\), *)?dtype=[^)]*\)",
    "= ...",
    stub,
  )

  stub = re.sub(r"\bpyvinecopulib\.[a-z_]+\.", "", stub)
  stub = re.sub(r"\bpyvinecopulib\.", "", stub)

  stub = re.sub(r"matplotlib\.figure\.Figure", "Figure", stub)
  stub = re.sub(r"matplotlib\.axes\._axes\.Axes", "Axes", stub)

  return stub


def generate_stub(
  module_name: str,
  output_path: Path,
  indent: int = 2,
) -> None:
  """Render a `.pyi` stub for ``module_name`` to ``output_path``.

  Assumes ``sys.path`` is already configured so that ``importlib.import_module``
  picks up the staged package.
  """
  pkg = importlib.import_module(module_name)
  names = sorted(getattr(pkg, "__all__", []))

  # Empty top-level __all__ ⇒ C++ extension not staged. Subpackages legitimately
  # have empty __all__ (e.g. sklearn placeholder).
  if not names and module_name == "pyvinecopulib":
    raise SystemExit(
      f"Refusing to overwrite {output_path}: imported {module_name} has an "
      f"empty __all__ (module path: {getattr(pkg, '__file__', '<namespace>')}).\n"
      "The C++ extension is likely not built into this environment. "
      "Run 'pip install -e . --no-build-isolation' and retry."
    )

  ext = importlib.import_module("pyvinecopulib.pyvinecopulib_ext")
  family_cls = ext.BicopFamily

  def _safe_getattr(obj_name):
    try:
      return getattr(pkg, obj_name, None)
    except ImportError:
      return None

  known_types = {name for name in names if inspect.isclass(_safe_getattr(name))}

  indent_str = " " * indent

  lines = [
    "import collections",
    "from typing import Any, Optional",
    "from numpy.typing import ArrayLike, NDArray",
    "from matplotlib.figure import Figure",
    "from matplotlib.axes import Axes",
    "",
  ]

  for name in names:
    # `getattr` may trigger a lazy submodule import that fails when the
    # corresponding extra isn't installed (e.g. sklearn on build machines).
    try:
      obj = getattr(pkg, name, None)
    except ImportError:
      lines.append(f"from . import {name} as {name}\n")
      continue

    if inspect.ismodule(obj):
      lines.append(f"from . import {name} as {name}\n")
      continue

    # Cross-subpackage canonical class: emit a relative re-export alias.
    canonical = obj.__module__ if inspect.isclass(obj) else None
    if (
      canonical
      and canonical != module_name
      and isinstance(canonical, str)
      and _is_same_package(canonical, module_name)
    ):
      relative_import = _relative_import_spec(canonical, module_name)
      if relative_import:
        lines.append(f"{relative_import} import {name} as {name}\n")
        continue

    if inspect.isclass(obj):
      lines.extend(
        render_class_stub(obj, name, known_types=known_types, indent=indent)
      )
    elif (
      isinstance(obj, (FunctionType, BuiltinFunctionType))
      or type(obj).__name__ == "nb_func"
    ):
      try:
        if (
          type(obj).__name__ == "nb_func"
          or hasattr(obj, "__module__")
          and obj.__module__
          and "pyvinecopulib_ext" in obj.__module__
        ):
          lines.extend(
            render_nanobind_function_stub(
              obj, name, known_types=known_types, indent=indent
            )
          )
        else:
          lines.extend(
            render_python_function_stub(
              obj, name, known_types=known_types, indent=indent
            )
          )
      except Exception as e:
        lines.append(
          f"{indent_str}# def {name}(...):  # signature unavailable ({e})"
        )
        lines.append(f"{indent_str}...\n")
    elif isinstance(obj, family_cls):
      lines.append(f"{name}: BicopFamily = ...\n")
    elif isinstance(obj, list) and all(isinstance(x, family_cls) for x in obj):
      lines.append(f"{name}: list[BicopFamily] = ...\n")
    elif name == "__version__" and isinstance(obj, str):
      lines.append("__version__: str = ...\n")
    else:
      lines.append(f"{name}: Any = ...\n")

  # Emit a permissive __getattr__ stub when the module has one, so static
  # checkers don't flag access to (e.g.) warn-on-access deprecated names.
  module_dict = vars(pkg)
  if "__getattr__" in module_dict and callable(module_dict["__getattr__"]):
    lines.append("def __getattr__(name: str) -> Any: ...\n")

  stub = cleanup_stub("\n".join(lines))
  output_path.parent.mkdir(parents=True, exist_ok=True)
  output_path.write_text(stub, encoding="utf-8")
  print(f"Wrote stub for {len(names)} symbols ({module_name}) to {output_path}")


def _is_same_package(target: str, current: str) -> bool:
  """Return True if target and current share a common top-level package."""
  return target.split(".", 1)[0] == current.split(".", 1)[0]


def _relative_import_spec(target: str, current: str) -> str:
  """Build a `from ...x.y` import prefix for `target` relative to `current`.

  Both inputs are dotted module paths. Returns "" if the import can't be
  expressed (e.g. fully unrelated packages).
  """
  t = target.split(".")
  c = current.split(".")
  # Find longest common prefix
  i = 0
  while i < len(t) and i < len(c) and t[i] == c[i]:
    i += 1
  if i == 0:
    return ""
  # Up `len(c) - i` levels, then down through t[i:]
  dots = "." * (len(c) - i + 1)
  tail = ".".join(t[i:])
  return f"from {dots}{tail}"


def _parse_module_output(spec: str) -> tuple[str, Path]:
  """Parse a ``PKG:PATH`` spec from --module-output."""
  if ":" not in spec:
    raise argparse.ArgumentTypeError(
      f"Expected --module-output in PKG:PATH form, got {spec!r}"
    )
  module_name, _, output = spec.partition(":")
  module_name = module_name.strip()
  output = output.strip()
  if not module_name or not output:
    raise argparse.ArgumentTypeError(
      f"--module-output requires non-empty PKG and PATH, got {spec!r}"
    )
  return module_name, Path(output)


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--source-pkg-dir",
    required=True,
    type=Path,
    help="Source pyvinecopulib/ directory (contains __init__.py).",
  )
  parser.add_argument(
    "--ext-so",
    required=True,
    type=Path,
    help="Path to the freshly built pyvinecopulib_ext shared library.",
  )
  parser.add_argument(
    "--module-output",
    action="append",
    default=[],
    type=_parse_module_output,
    metavar="PKG:PATH",
    help=(
      "Repeatable: emit a stub for module PKG at PATH. "
      "Modules whose import raises ImportError are skipped with a warning."
    ),
  )
  parser.add_argument(
    "--py-typed",
    type=Path,
    default=None,
    help="Optional path to py.typed marker (defaults to <source-pkg-dir>/py.typed).",
  )
  parser.add_argument(
    "--indent", type=int, default=2, help="Indentation level (spaces)."
  )
  args = parser.parse_args()

  if not args.module_output:
    sys.exit("At least one --module-output PKG:PATH is required.")

  if not args.ext_so.exists():
    sys.exit(f"Extension not found at: {args.ext_so}")
  if not args.source_pkg_dir.is_dir():
    sys.exit(f"Source package dir not found: {args.source_pkg_dir}")

  py_typed = args.py_typed or (args.source_pkg_dir / "py.typed")

  # Stage sources + freshly built .so/.pyd in a tempdir so importlib sees a
  # complete package. mkdtemp + best-effort rmtree (not TemporaryDirectory):
  # Windows holds the imported .pyd locked for the interpreter's lifetime,
  # which would make TemporaryDirectory's strict cleanup raise.
  tmp = tempfile.mkdtemp(prefix="pyvinecopulib-stubs-")
  try:
    site = Path(tmp)
    pkg = site / "pyvinecopulib"
    shutil.copytree(
      args.source_pkg_dir,
      pkg,
      ignore=shutil.ignore_patterns(
        "__init__.pyi", "py.typed", "*.so", "*.pyd", "*.dylib", "__pycache__"
      ),
    )
    shutil.copy2(args.ext_so, pkg / args.ext_so.name)

    sys.path.insert(0, str(site))
    try:
      for module_name, output_path in args.module_output:
        try:
          generate_stub(module_name, output_path, indent=args.indent)
        except ImportError as e:
          # Optional subpackages (e.g. sklearn without scikit-learn installed).
          print(
            f"Skipping {module_name}: import failed ({e}). "
            f"No stub written to {output_path}.",
            file=sys.stderr,
          )
    finally:
      sys.path.remove(str(site))

    py_typed.parent.mkdir(parents=True, exist_ok=True)
    py_typed.write_text("", encoding="utf-8")
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
  main()
