"""The writing rules codespell cannot express.

`.codespell-prose.txt` bans single tokens, which covers most of the rule and
runs in `make lint` and pre-commit. Two things are out of its reach and live
here instead:

* **Multi-word phrases.** codespell's word regex is `[\\w\\-'’]+`, so a phrase
  is several tokens and can never be one dictionary key.
* **The difference between prose and code.** codespell reads a file as text,
  so it cannot tell a docstring from an identifier -- which is why a `Tensor`
  named for a boolean mask has to be exempt, and why that exemption is pinned
  here rather than scattered as inline ignores.

A banned phrase is sometimes the right phrase. Wrap the lines in
``# codespell:ignore-begin`` / ``-end`` and both halves of the rule skip them,
which is the one escape hatch rather than two.
"""

from __future__ import annotations

import ast
import io
import pathlib
import re
import tokenize
from typing import Iterator

import pytest

#: Phrases that read as filler. Each is a regex over collapsed whitespace, so
#: a phrase split across two source lines is still found.
_BANNED_PHRASES = (
  r"the whole of",
  r"the whole point",
  r"the point is",
  r"keeps (?:\w+|it|that|the \w+) honest",
  r"worth not undoing",
  r"under the hood",
  r"at its core",
  r"the beauty of",
  r"deep dive",
  r"heavy lifting",
  r"in the wild",
  r"(?:precisely|exactly) why",
)

#: Files whose *identifiers* legitimately use a banned word, with the reason.
#: codespell has no way to skip an identifier while still reading the prose
#: around it, so the exemption is recorded rather than suppressed.
_IDENTIFIER_EXEMPTIONS = {
  "src/pyvinecopulib/torch/_fit_tll.py": "`gate` names a boolean mask, which "
  "is ordinary array-programming vocabulary; the numpydoc entry documents the "
  "parameter name, so renaming it would be a code change.",
}


def _sources() -> Iterator[pathlib.Path]:
  """Every Python file whose prose this rule covers.

  Yields
  ------
  pathlib.Path
      A source file under ``src`` or ``tests``.
  """
  for root in ("src", "tests"):
    base = pathlib.Path(root)
    if base.is_dir():
      for path in sorted(base.rglob("*.py")):
        # This file has to spell the words it bans, so it cannot be its own
        # subject -- the same reason the AGENTS.md rule points at the
        # dictionary instead of quoting it.
        if path.name != "test_prose.py":
          yield path


def _without_exempt_regions(text: str) -> str:
  """Blank out the lines inside a ``codespell:ignore`` region.

  The same escape hatch the single-token half of the rule uses, honored here so
  a phrase and a word are exempted the same way. Lines are blanked rather than
  removed so every reported line number still points at the real line.

  Parameters
  ----------
  text : str
      The file's source.

  Returns
  -------
  str
      The source with exempt regions replaced by blank lines.
  """
  out: list[str] = []
  inside = False
  for line in text.splitlines(True):
    if "codespell:ignore-begin" in line:
      inside = True
    elif "codespell:ignore-end" in line:
      inside = False
      out.append("\n")
      continue
    out.append("\n" if inside else line)
  return "".join(out)


def _prose_of(path: pathlib.Path) -> list[tuple[int, str]]:
  """The docstrings and comments in one file, as ``(line, text)``.

  Parameters
  ----------
  path : pathlib.Path
      The file to read.

  Returns
  -------
  list of tuple
      One entry per docstring or comment.
  """
  src = _without_exempt_regions(path.read_text(encoding="utf-8"))
  out: list[tuple[int, str]] = []
  tree = ast.parse(src)
  for node in ast.walk(tree):
    if isinstance(
      node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
    ):
      doc = ast.get_docstring(node, clean=False)
      if doc:
        out.append((getattr(node, "lineno", 1), doc))
  for tok in tokenize.generate_tokens(io.StringIO(src).readline):
    if tok.type == tokenize.COMMENT:
      out.append((tok.start[0], tok.string))
  return out


@pytest.mark.parametrize("phrase", _BANNED_PHRASES)
def test_no_banned_phrase_in_prose(phrase: str) -> None:
  """A phrase codespell cannot tokenize is still banned.

  Whitespace is collapsed first, so wrapping a phrase across two lines does
  not hide it -- which is how several of these were written in the first
  place.
  """
  pattern = re.compile(phrase, re.IGNORECASE)
  found: list[str] = []
  for path in _sources():
    for line, text in _prose_of(path):
      if pattern.search(" ".join(text.split())):
        found.append(f"{path}:{line}")
  assert found == [], f"{phrase!r} appears at: {found}"


def test_the_identifier_exemptions_are_still_needed() -> None:
  """An exemption outlives its reason unless something checks.

  Each entry claims a file uses a banned word as an *identifier*. If the
  identifier is gone, the entry should be too -- otherwise the list slowly
  becomes a place to hide new violations.
  """
  for rel, reason in _IDENTIFIER_EXEMPTIONS.items():
    path = pathlib.Path(rel)
    if not path.is_file():  # installed rather than checked out
      pytest.skip(f"{rel} not present")
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = {
      node.arg for node in ast.walk(tree) if isinstance(node, ast.arg)
    } | {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
    assert names & {"gate"}, f"{rel}: exemption no longer applies -- {reason}"


def test_the_ban_list_is_wired_into_codespell() -> None:
  """The single-token half of the rule has to actually run.

  Two mechanics are easy to get wrong and silent when wrong: the dictionary
  file must come *first* in the setting, because a value beginning with ``-``
  is read as another flag, and the trailing ``-`` is what keeps the
  ``en-GB_to_en-US`` builtins the American-English rule depends on.
  """
  pyproject = pathlib.Path("pyproject.toml")
  banlist = pathlib.Path(".codespell-prose.txt")
  if not pyproject.is_file():  # installed rather than checked out
    pytest.skip("source tree not available")

  text = pyproject.read_text(encoding="utf-8")
  assert 'dictionary = ".codespell-prose.txt,-"' in text
  assert 'builtin = "clear,rare,en-GB_to_en-US"' in text
  assert banlist.is_file()

  # Every entry must carry a reason after a comma: that is what makes
  # codespell report a suggestion without rewriting the word, since the right
  # replacement depends on the sentence.
  for raw in banlist.read_text(encoding="utf-8").splitlines():
    if not raw.strip():
      continue
    assert "->" in raw, raw
    _, data = raw.split("->", 1)
    assert "," in data, f"no reason, so codespell would auto-fix it: {raw!r}"
