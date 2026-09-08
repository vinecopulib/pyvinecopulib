"""The worked examples in docstrings are executed, not just written.

AGENTS.md requires every code example to run somewhere: a notebook cell or a
doctest. The four contract examples in ``core/protocols.py`` are neither --
they are ``::`` literal blocks appended to both a protocol's docstring and its
canonical base's, so they render on the public docs for eight classes while
nothing ran them. `VinedistLike`'s opened with ``to_pseudo_obs(y)`` and never
defined ``y``.

Doctest would be the idiomatic home and is not wired up; until it is, this is
what keeps the promise.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import re
import textwrap

import pytest

#: The module whose ``_*_EXAMPLE`` blocks are rendered on the docs site.
_SOURCE = pathlib.Path("src/pyvinecopulib/core/protocols.py")


def _literal_blocks(body: str) -> list[str]:
  """The reStructuredText literal blocks introduced by a ``::`` line.

  A block runs to the first line indented less than its own first line, which
  is what separates it from the prose that follows -- several of these examples
  have more prose after the code.

  Parameters
  ----------
  body : str
      The docstring fragment to scan.

  Returns
  -------
  list of str
      One dedented, executable block per ``::``.
  """
  out: list[str] = []
  lines = body.splitlines(True)
  i = 0
  while i < len(lines):
    if not lines[i].rstrip().endswith("::"):
      i += 1
      continue
    j = i + 1
    while j < len(lines) and not lines[j].strip():
      j += 1
    if j >= len(lines):
      break
    indent = len(lines[j]) - len(lines[j].lstrip())
    chunk: list[str] = []
    while j < len(lines):
      line = lines[j]
      if not line.strip():
        chunk.append("\n")
      elif len(line) - len(line.lstrip()) >= indent:
        chunk.append(line)
      else:
        break
      j += 1
    out.append(textwrap.dedent("".join(chunk)))
    i = j
  return out


def _examples() -> dict[str, list[str]]:
  """Every ``_*_EXAMPLE`` in the source, as executable blocks.

  Returns
  -------
  dict
      Name to its list of blocks.
  """
  src = _SOURCE.read_text(encoding="utf-8")
  found: dict[str, list[str]] = {}
  for name in re.findall(r"^(_\w+_EXAMPLE) = ", src, re.M):
    match = re.search(rf'^{name} = """(.*?)^"""', src, re.S | re.M)
    assert match is not None, name
    found[name] = _literal_blocks(match.group(1))
  return found


@pytest.mark.parametrize(
  "name",
  [
    "_BICOP_EXAMPLE",
    "_VINECOP_EXAMPLE",
    "_MARGIN_EXAMPLE",
    "_VINEDIST_EXAMPLE",
  ],
)
def test_a_contracts_worked_example_runs(name: str) -> None:
  """Run one contract's example, blocks sharing a namespace in order.

  The later blocks of an example continue the earlier ones -- `VinedistLike`'s
  second block refits the distribution built by its first -- so they run in
  one namespace rather than independently.
  """
  if not _SOURCE.is_file():  # installed rather than checked out
    pytest.skip("source tree not available")
  blocks = _examples()[name]
  assert blocks, f"{name} has no literal block"
  if any("scipy" in block for block in blocks):
    pytest.importorskip("scipy.stats")

  namespace: dict[str, object] = {}
  for index, block in enumerate(blocks):
    with contextlib.redirect_stdout(io.StringIO()):
      exec(compile(block, f"{name}#{index}", "exec"), namespace)


def test_every_example_in_the_module_is_covered() -> None:
  """The parametrization is a list, so it has to be checked against the source.

  A fifth contract would otherwise render an unexecuted example, which is the
  state this file exists to end.
  """
  if not _SOURCE.is_file():  # installed rather than checked out
    pytest.skip("source tree not available")
  covered = {
    "_BICOP_EXAMPLE",
    "_VINECOP_EXAMPLE",
    "_MARGIN_EXAMPLE",
    "_VINEDIST_EXAMPLE",
  }
  assert set(_examples()) == covered
