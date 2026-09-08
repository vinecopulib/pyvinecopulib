"""The writing rules codespell cannot express.

`.codespell-prose.txt` bans single tokens, which covers most of the rule and
runs in `make lint` and pre-commit. Two things are out of its reach and live
here instead:

* **Multi-word phrases.** codespell's word regex is `[\\w\\-'’]+`, so a phrase
  is several tokens and can never be one dictionary key.
* **Banned words inside identifiers.** That same regex makes
  `_GENUINE_OPT_OUTS` a *single* token, which matches no dictionary key, so
  the whole name is invisible to codespell. This file reads the identifiers
  out of the AST and checks them against the same list.
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
  # Multi-word, so codespell's tokenizer can never hold them.
  r"byte for byte",
  r"escape hatch",
  r"blind sweep",
  r"in flight",
)

#: Files whose *identifiers* legitimately use a banned word, with the reason.
#: codespell has no way to skip an identifier while still reading the prose
#: around it, so the exemption is recorded rather than suppressed.
_IDENTIFIER_EXEMPTIONS: dict[str, str] = {}


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


def _banned_words() -> set[str]:
  """The single-token half of the rule, read from the dictionary itself.

  Returns
  -------
  set of str
      Every banned word, lowercased.
  """
  path = pathlib.Path(".codespell-prose.txt")
  if not path.is_file():  # installed rather than checked out
    return set()
  out: set[str] = set()
  for raw in path.read_text(encoding="utf-8").splitlines():
    if "->" in raw:
      out.add(raw.split("->", 1)[0].strip().lower())
  return out


def _identifiers_of(path: pathlib.Path) -> Iterator[tuple[int, str]]:
  """Every name a file binds or reads, as ``(line, name)``.

  Yields
  ------
  tuple
      The line and the identifier.
  """
  tree = ast.parse(_without_exempt_regions(path.read_text(encoding="utf-8")))
  for node in ast.walk(tree):
    if isinstance(node, ast.Name):
      yield node.lineno, node.id
    elif isinstance(node, ast.arg):
      yield node.lineno, node.arg
    elif isinstance(
      node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
    ):
      yield node.lineno, node.name
    elif isinstance(node, ast.Attribute):
      yield node.lineno, node.attr


def test_no_banned_word_hides_inside_an_identifier() -> None:
  """`_GENUINE_OPT_OUTS` passed codespell because it is one token.

  The dictionary is the source of truth, so this reads it rather than
  restating it. An identifier that legitimately uses a banned word is
  recorded in ``_IDENTIFIER_EXEMPTIONS`` or wrapped in an ignore region --
  the same two escapes the rest of the rule offers.
  """
  banned = _banned_words()
  if not banned:
    pytest.skip("ban list not present")
  exempt = {pathlib.Path(rel) for rel in _IDENTIFIER_EXEMPTIONS}

  found: list[str] = []
  for path in _sources():
    if path in exempt:
      continue
    for line, name in _identifiers_of(path):
      parts = {p.lower() for p in re.split(r"[^A-Za-z]+", name) if p}
      hit = parts & banned
      if hit:
        found.append(f"{path}:{line}: {name} contains {sorted(hit)}")
  assert found == [], found


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

  Each entry claims a file needs a banned word as an *identifier*. The list is
  empty: every case so far had a name that read better anyway -- the mask
  `_fit_tll.py` called `gate` is the `outer` condition its own docstring
  describes. Keep it that way if you can; the entry is the fallback.
  """
  for rel, reason in _IDENTIFIER_EXEMPTIONS.items():
    path = pathlib.Path(rel)
    if not path.is_file():  # installed rather than checked out
      pytest.skip(f"{rel} not present")
    banned = _banned_words()
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = {
      node.arg for node in ast.walk(tree) if isinstance(node, ast.arg)
    } | {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
    parts = {
      part.lower() for name in names for part in re.split(r"[^A-Za-z]+", name)
    }
    assert parts & banned, f"{rel}: exemption no longer applies -- {reason}"


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
