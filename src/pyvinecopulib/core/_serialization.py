"""JSON persistence for the margin and distribution layer.

``Bicop`` / ``Vinecop`` / ``RVineStructure`` carry ``to_json`` / ``from_json``
from the compiled library, and ``docs/concepts.rst`` presents that as the way to
persist a model. The objects this package adds on top -- ``Vinedist`` and the
margins -- need the same surface, and a margin can be any class a caller wrote,
so reconstruction goes through a registry rather than a fixed list of types.

The payload carries a ``kind`` naming the margin's type and a ``version`` that a
reader checks, so a format change is a loud failure rather than a wrong model.
"""

from __future__ import annotations

import json
import math
from typing import Any, Callable

#: Bumped when the payload's shape changes incompatibly.
MARGIN_JSON_VERSION = 1

#: ``kind`` -> a callable rebuilding the margin from its payload.
_READERS: dict[str, Callable[[dict[str, Any]], Any]] = {}


def register_margin_json(
  kind: str, reader: Callable[[dict[str, Any]], Any]
) -> None:
  """Teach :func:`margin_from_json` how to rebuild one margin type.

  Parameters
  ----------
  kind : str
      The value the margin's ``to_json`` writes under ``"kind"``. Conventionally
      the class name.
  reader : callable
      ``reader(payload) -> margin``, receiving the mapping that
      :func:`margin_to_json` produced.

  Raises
  ------
  ValueError
      If ``kind`` is already registered to a different reader.
  """
  existing = _READERS.get(kind)
  if existing is not None and existing is not reader:
    raise ValueError(f"a reader for margin kind {kind!r} is already registered")
  _READERS[kind] = reader


def margin_to_json(margin: Any) -> dict[str, Any]:
  """Return one margin's JSON payload.

  Parameters
  ----------
  margin : MarginLike
      The margin to serialize. It must provide ``to_json``, which every margin
      this package ships does.

  Returns
  -------
  dict
      A JSON-serializable mapping carrying ``kind`` and ``version``.

  Raises
  ------
  TypeError
      If the margin has no ``to_json``.
  """
  to_json = getattr(margin, "to_json", None)
  if to_json is None:
    raise TypeError(
      f"{type(margin).__name__} cannot be serialized: it has no `to_json`. "
      "Implement `to_json` returning a JSON-serializable mapping, and call "
      "`pyvinecopulib.core.register_margin_json` so it can be read back."
    )
  payload = to_json()
  if isinstance(payload, str):
    # A compiled margin (`Kde1d`) serializes itself to a JSON string.
    payload = {"kind": type(margin).__name__, "json": payload}
  else:
    payload = dict(payload)
    payload.setdefault("kind", type(margin).__name__)
  payload.setdefault("version", MARGIN_JSON_VERSION)
  return payload


def margin_from_json(payload: dict[str, Any]) -> Any:
  """Rebuild one margin from the payload :func:`margin_to_json` produced.

  Parameters
  ----------
  payload : dict
      A mapping carrying ``kind`` and ``version``.

  Returns
  -------
  MarginLike
      The reconstructed margin.

  Raises
  ------
  ValueError
      If ``kind`` is unknown or unregistered, or the version is unrecognized.
  """
  kind = payload.get("kind")
  version = payload.get("version")
  if version != MARGIN_JSON_VERSION:
    raise ValueError(
      f"unsupported margin JSON version {version!r}; this build reads "
      f"version {MARGIN_JSON_VERSION}"
    )
  reader = _READERS.get(kind)
  if reader is None:
    known = ", ".join(sorted(_READERS)) or "(none)"
    raise ValueError(
      f"no reader registered for margin kind {kind!r}; known kinds: {known}. "
      "Call `pyvinecopulib.core.register_margin_json` first."
    )
  return reader(payload)


def _encode_nonfinite(value: Any) -> Any:
  """Replace non-finite floats with strings, recursively.

  ``json.dumps`` writes ``Infinity`` / ``NaN``, which strict JSON has no
  literal for -- and which the C++ reader behind :func:`write_file` rejects.
  They travel as strings and are restored on read, so a ``-inf`` log-likelihood
  in a selection report survives exactly rather than becoming ``null``.

  Parameters
  ----------
  value : object
      Any JSON-serializable structure.

  Returns
  -------
  object
      The same structure with non-finite floats replaced.
  """
  if isinstance(value, float) and not math.isfinite(value):
    return f"__nonfinite__:{value!r}"
  if isinstance(value, dict):
    return {k: _encode_nonfinite(v) for k, v in value.items()}
  if isinstance(value, (list, tuple)):
    return [_encode_nonfinite(v) for v in value]
  return value


def _decode_nonfinite(value: Any) -> Any:
  """Invert :func:`_encode_nonfinite`.

  Parameters
  ----------
  value : object
      A structure parsed from JSON.

  Returns
  -------
  object
      The same structure with the encoded floats restored.
  """
  if isinstance(value, str) and value.startswith("__nonfinite__:"):
    return float(value.split(":", 1)[1])
  if isinstance(value, dict):
    return {k: _decode_nonfinite(v) for k, v in value.items()}
  if isinstance(value, list):
    return [_decode_nonfinite(v) for v in value]
  return value


def dumps(payload: dict[str, Any]) -> str:
  """Serialize a payload to a JSON string.

  Parameters
  ----------
  payload : dict
      The mapping to serialize.

  Returns
  -------
  str
      Its JSON representation, with non-finite floats encoded as strings so the
      result is strict JSON.
  """
  return json.dumps(_encode_nonfinite(payload), allow_nan=False)


def loads(text: str) -> dict[str, Any]:
  """Parse a JSON string into a payload.

  Parameters
  ----------
  text : str
      A JSON object produced by :func:`dumps`.

  Returns
  -------
  dict
      The parsed mapping.

  Raises
  ------
  ValueError
      If the text is not a JSON object.
  """
  payload = _decode_nonfinite(json.loads(text))
  if not isinstance(payload, dict):
    raise ValueError("expected a JSON object")
  return payload


def _register_builtin_readers() -> None:
  """Register the readers for the margins this package ships."""

  def _kde1d(payload: dict[str, Any]) -> Any:
    from . import Kde1d

    return Kde1d.from_json(payload["json"])

  def _parametric(payload: dict[str, Any]) -> Any:
    from ..margins import ParametricMargin

    return ParametricMargin._from_json_payload(payload)

  def _selector(payload: dict[str, Any]) -> Any:
    from ..margins import MarginSelector

    return MarginSelector._from_json_payload(payload)

  register_margin_json("Kde1d", _kde1d)
  register_margin_json("ParametricMargin", _parametric)
  register_margin_json("MarginSelector", _selector)


_register_builtin_readers()


def write_file(filename: str, text: str) -> None:
  """Write a JSON payload, as CBOR when the name ends in ``.cbor``.

  The extension rule is the one ``Bicop.to_file`` / ``Vinecop.to_file`` follow
  -- the same C++ helper, so the whole model surface reads and writes the same
  formats.

  Parameters
  ----------
  filename : str
      Path to write.
  text : str
      A JSON string.
  """
  from ..pyvinecopulib_ext import _json_to_file

  _json_to_file(filename, text)


def read_file(filename: str) -> str:
  """Read a JSON payload written by :func:`write_file`.

  Parameters
  ----------
  filename : str
      Path to read.

  Returns
  -------
  str
      The payload as a JSON string.
  """
  from ..pyvinecopulib_ext import _file_to_json

  return _file_to_json(filename)
