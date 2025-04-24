#!/usr/bin/env python3
"""
Generate __init__.pyi from __all__
"""

import argparse
import importlib
import inspect
import re
import sys
from pathlib import Path
from types import BuiltinFunctionType, FunctionType
from typing import Optional


def wrap_known_types(sig: str, known_types: set[str]) -> str:
  for name in known_types:
    # Replace list of constants like [BicopFamily.indep, BicopFamily.gaussian, ...]
    sig = re.sub(
      rf"= *\[{name}\.\w+(?: *, *{name}\.\w+)*\]",
      "= ...",
      sig,
    )

    # Replace default constructor calls like FitControlsBicop() with ...
    sig = re.sub(rf"= *{name}\s*\(\)", "= ...", sig)

    # Replace constant access like BicopFamily.indep with ...
    sig = re.sub(rf"= *{name}\s*\.\w+", "= ...", sig)

    # Wrap standalone type names (not followed by '(' or '.')
    sig = re.sub(rf"\b{name}\b(?!\s*\()", f'"{name}"', sig)

    # # # Wrap standalone type names (except when followed by a '.')
    # sig = re.sub(rf"\b{name}\b(?!\s*\.)", f'"{name}"', sig)

  return sig


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
    except Exception:
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
  # Replace any annotation like numpy.ndarray[...] with ArrayLike
  stub = re.sub(r"numpy\.ndarray\[.*?\]", "ArrayLike", stub)

  # Replace default values like = array([], dtype=float64) or np.array([], dtype=np.float64)
  stub = re.sub(r"= *(?:np\.)?array\(\[\], *dtype=.*?\)", "= ...", stub)
  stub = re.sub(
    r"= *(?:np\.)?array\(\[\], *(?:shape=\([^\)]*\), *)?dtype=[^)]*\)",
    "= ...",
    stub,
  )

  # # Replace FitControlsBicop() and FitControlsVinecop() with ...
  # stub = re.sub(r"= *FitControlsBicop\(\)", "= ...", stub)
  # stub = re.sub(r"= *FitControlsVinecop\(\)", "= ...", stub)

  # Remove pyvinecopulib. prefix from types
  stub = re.sub(r"\bpyvinecopulib\.", "", stub)

  return stub


def generate_stub(site_dir: str, output_path: Path, indent: int = 2):
  sys.path.insert(0, site_dir)
  pkg = importlib.import_module("pyvinecopulib")
  names = sorted(getattr(pkg, "__all__", []))

  known_types = {
    name for name in names if inspect.isclass(getattr(pkg, name, None))
  }

  indent_str = " " * indent

  lines = [
    "import collections",
    "from typing import Any",
    "from numpy.typing import ArrayLike",
    "",
  ]

  for name in names:
    obj = getattr(pkg, name, None)

    if inspect.isclass(obj):
      lines.extend(
        render_class_stub(obj, name, known_types=known_types, indent=indent)
      )
    elif (
      isinstance(obj, (FunctionType, BuiltinFunctionType))
      or type(obj).__name__ == "nb_func"
    ):
      try:
        lines.extend(
          render_nanobind_function_stub(
            obj, name, known_types=known_types, indent=indent
          )
        )
      except Exception as e:
        lines.append(
          f"{indent_str}# def {name}(...):  # signature unavailable ({e})"
        )
        lines.append(f"{indent_str}...\n")
    elif isinstance(obj, pkg.BicopFamily):
      lines.append(f"{name}: BicopFamily = ...\n")
    elif isinstance(obj, list) and all(
      isinstance(x, pkg.BicopFamily) for x in obj
    ):
      lines.append(f"{name}: list[BicopFamily] = ...\n")
    elif name == "__version__" and isinstance(obj, str):
      lines.append("__version__: str = ...\n")
    else:
      lines.append(f"{name}: Any = ...\n")

  stub = cleanup_stub("\n".join(lines))
  output_path.write_text(stub, encoding="utf-8")
  print(f"Wrote stub for {len(names)} symbols to {output_path}")


def main():
  parser = argparse.ArgumentParser(
    description="Generate .pyi stub for pyvinecopulib from __all__."
  )
  parser.add_argument("site_dir", help="Path to site-packages directory")
  parser.add_argument(
    "--indent", type=int, default=2, help="Indentation level (spaces)"
  )
  args = parser.parse_args()

  output_path = Path("src/pyvinecopulib/__init__.pyi")
  generate_stub(args.site_dir, output_path, indent=args.indent)
  Path("src/pyvinecopulib/py.typed").write_text("", encoding="utf-8")


if __name__ == "__main__":
  main()
