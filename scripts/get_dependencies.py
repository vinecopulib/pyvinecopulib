import sys
from pathlib import Path

if sys.version_info >= (3, 11):
  import tomllib
else:
  try:
    import tomli as tomllib
  except ImportError:
    print(
      "tomli is required for Python < 3.11. Run: pip install tomli",
      file=sys.stderr,
    )
    sys.exit(1)


def parse_pyproject_toml() -> dict:
  pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
  return tomllib.loads(pyproject_path.read_text())


def get_dependencies(
  groups: list[str], include_core: bool = False, include_build: bool = False
) -> list[str]:
  data = parse_pyproject_toml()

  project = data.get("project", {})
  build_sys = data.get("build-system", {})

  extras = project.get("optional-dependencies", {})
  core = project.get("dependencies", [])
  build = build_sys.get("requires", [])

  result = []

  if include_core:
    result.extend(core)

  if include_build:
    result.extend(build)

  for group in groups:
    if group not in extras:
      raise ValueError(f"Dependency group '{group}' not found.")
    result.extend(extras[group])

  return result


def main():
  import argparse

  parser = argparse.ArgumentParser(
    description="Extract dependencies from pyproject.toml"
  )
  parser.add_argument(
    "groups", nargs="*", help="Optional dependency groups to include"
  )
  parser.add_argument(
    "--include-core", action="store_true", help="Include core dependencies"
  )
  parser.add_argument(
    "--include-build",
    action="store_true",
    help="Include build dependencies from [build-system.requires]",
  )
  parser.add_argument(
    "--as-requirements",
    action="store_true",
    help="Output dependencies one per line (requirements.txt style)",
  )
  args = parser.parse_args()

  try:
    deps = get_dependencies(
      args.groups,
      include_core=args.include_core,
      include_build=args.include_build,
    )
  except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)

  if args.as_requirements:
    for dep in deps:
      print(dep.strip())
  else:
    print(" ".join(dep.strip() for dep in deps))


if __name__ == "__main__":
  main()
