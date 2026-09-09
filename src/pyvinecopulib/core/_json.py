"""How a model payload is encoded and written.

Not about margins: :meth:`~pyvinecopulib.core.VinedistBase.to_json` stamps a
whole distribution -- its copula and every margin -- through the same
functions, so they live beside neither half.

The payload carries a ``version`` a reader checks, and non-finite floats
travel as strings because JSON has no spelling for them.
"""

from __future__ import annotations

import json
import math
from typing import Any, cast

__all__ = ["MODEL_JSON_VERSION", "dumps", "loads"]

#: Bumped when a payload's shape changes incompatibly. One number for the
#: whole model: a margin payload and the distribution payload that embeds it
#: are read by the same build, so a change to either breaks both.
MODEL_JSON_VERSION = 1


def _encode_nonfinite(value: object) -> object:
  """Replace non-finite floats with strings, recursively.

  ``json.dumps`` writes ``Infinity`` / ``NaN``, which strict JSON has no
  literal for -- and which the reader behind :func:`write_file` rejects.
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


def _decode_nonfinite(value: object) -> object:
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
      A JSON object produced by ``dumps``.

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


def write_file(filename: str, text: str) -> None:
  """Write a JSON payload, as CBOR when the name ends in ``.cbor``.

  The extension rule is the one ``Bicop.to_file`` / ``Vinecop.to_file`` follow
  -- the same helper, so the whole model surface reads and writes the same
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

  return cast("str", _file_to_json(filename))
