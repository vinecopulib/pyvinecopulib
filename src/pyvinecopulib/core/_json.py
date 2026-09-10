"""How a model payload is encoded and written.

Not about margins: :meth:`~pyvinecopulib.core.VinedistBase.to_json` stamps a
whole distribution -- its copula and every margin -- through the same
functions, so they live beside neither half.

The payload carries a ``version`` a reader checks, and a non-finite float
travels as a one-key tagged object because JSON has no spelling for one.
"""

from __future__ import annotations

import json
import math
from typing import Any, Optional, Union, cast

__all__ = ["MODEL_JSON_VERSION", "dumps", "loads", "read_payload"]

#: Bumped when a payload's shape changes incompatibly. One number for the
#: whole model: a margin payload and the distribution payload that embeds it
#: are read by the same build, so a change to either breaks both.
MODEL_JSON_VERSION = 1


#: Single key marking a tagged non-finite float, namespaced so no payload
#: field collides with it.
_NONFINITE_KEY = "__pyvinecopulib_nonfinite__"

#: The closed set of tags, which is what keeps the decode from being handed
#: arbitrary text: anything else stays the mapping it arrived as.
_NONFINITE_VALUES = {"nan": math.nan, "+inf": math.inf, "-inf": -math.inf}


def _encode_nonfinite(value: object) -> object:
  """Replace non-finite floats with a tagged object, recursively.

  ``json.dumps`` writes ``Infinity`` / ``NaN``, which strict JSON has no
  literal for -- and which the reader behind :func:`write_file` rejects. So
  one travels as ``{_NONFINITE_KEY: tag}`` and is restored on read, and a
  ``-inf`` log-likelihood in a selection report survives exactly rather than
  becoming ``null``.

  A one-key object rather than a marked *string*, because payloads carry
  arbitrary user text -- a variable's name, a margin's ``family_name`` -- and
  a marked string is a value user data can accidentally spell. One that did
  came back a float.

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
    if math.isnan(value):
      tag = "nan"
    else:
      tag = "+inf" if value > 0 else "-inf"
    return {_NONFINITE_KEY: tag}
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
  if (
    isinstance(value, dict)
    and len(value) == 1
    and value.get(_NONFINITE_KEY) in _NONFINITE_VALUES
  ):
    return _NONFINITE_VALUES[cast("str", value[_NONFINITE_KEY])]
  # An unrecognized tag stays a mapping rather than raising: a payload written
  # by a newer build is caught by the version check in `read_payload`, which
  # can say so, and this cannot.
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


def read_payload(
  payload: Union[str, dict[str, Any]],
  what: str,
  *,
  kind: Optional[str] = None,
) -> dict[str, Any]:
  """Decode a payload if it is still a string, and check what it claims to be.

  The three checks every reader owes, in one place: a payload is decoded
  through :func:`loads` rather than ``json.loads`` so a non-finite float comes
  back as one, its ``version`` is compared against this build's, and its
  ``kind`` -- where the caller knows which class it is reading -- against the
  class asked for. A subclass's payload read as its base is quietly the wrong
  model, since a subclass is a different distribution.

  Parameters
  ----------
  payload : str or dict
      A JSON string, or an already-decoded mapping.
  what : str
      Names the thing being read, for the version message.
  kind : str, or None, optional
      Class name the payload must claim. ``None`` skips the check, for a
      reader that dispatches on ``kind`` itself and can say more about it.

  Returns
  -------
  dict
      The decoded payload.

  Raises
  ------
  ValueError
      If the version is unrecognized, or ``kind`` names another class.
  """
  decoded = loads(payload) if isinstance(payload, str) else dict(payload)
  version = decoded.get("version")
  if version != MODEL_JSON_VERSION:
    raise ValueError(
      f"unsupported {what} JSON version {version!r}; this build reads "
      f"version {MODEL_JSON_VERSION}"
    )
  claimed = decoded.get("kind")
  if kind is not None and claimed != kind:
    raise ValueError(
      f"this payload was written by {claimed!r}, not {kind!r}; read it back "
      f"with {claimed}.from_json, or write it with {kind}.to_json"
    )
  return decoded
