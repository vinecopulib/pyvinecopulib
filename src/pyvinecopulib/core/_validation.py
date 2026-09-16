"""The package's validators, so one rule and one message serve every layer.

Every check a margin, pair copula, vine or vine distribution performs on its
inputs lives here: the univariate and covariate layouts, observation weights,
and the fit-time refusal of covariates a part cannot read. Each takes ``name=``
so the message names the argument the caller actually passed. The refusal an
absent optional dependency earns lives here too, so every extra is named the
same way.
"""

from __future__ import annotations

import contextlib
from typing import Any, Iterator, Optional, cast

from .protocols import array_namespace

from .protocols import ArrayT

__all__ = [
  "reject_covariates",
  "usable_observations",
  "validate_weights",
]


def check_var_types(var_types: Optional[list[str]], d: int) -> tuple[str, ...]:
  """Normalize and validate a vine's per-variable types.

  Parameters
  ----------
  var_types : list of str, or None, optional
      Per-variable types, ``"c"`` or ``"d"``; ``None`` means all continuous.
  d : int
      Dimension the types must cover.

  Returns
  -------
  tuple of str
      The validated types, one per variable.

  Raises
  ------
  ValueError
      If the length is not ``d`` or an entry is outside ``{"c", "d"}``.
  """
  types = ("c",) * d if var_types is None else tuple(var_types)
  if len(types) != d:
    raise ValueError(f"var_types has {len(types)} entries, expected {d}")
  bad = [t for t in types if t not in ("c", "d")]
  if bad:
    raise ValueError(f"var_types entries must be 'c' or 'd'; got {bad[0]!r}")
  return types


def validate_univariate(values: ArrayT, *, name: str = "y") -> ArrayT:
  """Require the documented one-dimensional univariate-data layout.

  Parameters
  ----------
  values : array
      The data to check.
  name : str, default="y"
      The argument's name, used in the error message.

  Returns
  -------
  array
      ``values`` unchanged.

  Raises
  ------
  ValueError
      If ``values`` is not one-dimensional.
  """
  if getattr(values, "ndim", None) != 1:
    shape = tuple(getattr(values, "shape", ()))
    raise ValueError(f"{name} must have shape (n,); got {shape}")
  return values


def validate_covariates(
  x: Optional[ArrayT], n_rows: int, *, name: str = "x"
) -> None:
  """Require a two-dimensional covariate matrix row-aligned with the data.

  Parameters
  ----------
  x : array, or None, optional
      The covariate matrix, or ``None`` to accept.
  n_rows : int
      Number of observations the covariates must align with.
  name : str, default="x"
      The argument's name, used in the error message.

  Returns
  -------
  None

  Raises
  ------
  ValueError
      If ``x`` is not two-dimensional or has the wrong number of rows.
  """
  if x is None:
    return
  if getattr(x, "ndim", None) != 2:
    shape = tuple(getattr(x, "shape", ()))
    one_d = getattr(x, "ndim", None) == 1
    raise ValueError(
      f"{name} must have shape (n, p), with one row per observation; "
      f"got {shape}"
      + (
        ". A one-dimensional x is refused because it says nothing about "
        "which axis is which: reshape it to (n, 1) for one covariate per "
        "observation, or use `covariate_row(x)` for one row shared across "
        "them."
        if one_d
        else ""
      )
    )
  rows = int(cast("Any", x).shape[0])
  if rows != n_rows:
    raise ValueError(
      f"{name} must have one row per observation; got {rows} rows for "
      f"{n_rows} observations"
    )


def validate_weights(
  weights: Optional[ArrayT], values: ArrayT, *, name: str = "weights"
) -> Optional[ArrayT]:
  """Normalize and validate one real, finite, nonnegative weight per row.

  Parameters
  ----------
  weights : array, or None, optional
      The weights to check, or ``None`` to accept.
  values : array
      The observations the weights align with; also fixes the array namespace
      and device the weights are coerced onto.
  name : str, default="weights"
      The argument's name, used in the error message.

  Returns
  -------
  array, or None
      The weights on the observations' namespace, or ``None``.

  Notes
  -----
  ``NaN`` and ``0`` mark a *dropped* observation, the same convention the data
  follow -- so a weight vector may carry either, and what must hold is that
  something survives them. ``+/-inf`` and a negative weight are not drop
  markers and are refused.

  Raises
  ------
  ValueError
      If the weights are not one per observation, are infinite or negative, or
      leave no observation standing.
  TypeError
      If they do not have a real numeric dtype.
  """
  if weights is None:
    return None
  xp = array_namespace(values)
  try:
    weights = xp.asarray(
      weights,
      device=getattr(values, "device", None),
    )
  except (TypeError, ValueError) as exc:
    raise ValueError(
      f"{name} must be an array compatible with the observations"
    ) from exc
  ref: Any = values
  if weights.ndim != 1 or int(weights.shape[0]) != int(ref.shape[0]):
    raise ValueError(
      f"{name} must have shape ({int(ref.shape[0])},), with one weight "
      f"per observation; got {tuple(weights.shape)}"
    )
  if not xp.isdtype(weights.dtype, ("real floating", "integral")):
    raise TypeError(
      f"{name} must have a real numeric dtype; got {weights.dtype}"
    )
  if bool(xp.any(xp.isinf(weights))):
    raise ValueError(f"{name} must not contain infinite values")
  # `NaN < 0` is False, so this reads only the entries that are not drops.
  if bool(xp.any(weights < 0)):
    raise ValueError(f"{name} must be nonnegative")
  # Every entry being a drop marker leaves nothing to fit, which `Kde1d` and
  # `utils.wdm` both refuse. Checked here too so that the message names the
  # argument the caller passed, and so that a margin that delegates to neither
  # is held to the same rule.
  kept = weights[~xp.isnan(weights)]
  if not float(xp.sum(kept)) > 0.0:
    raise ValueError(
      f"{name} must leave at least one observation standing; NaN and 0 mark a "
      "dropped observation, and every entry is one"
    )
  return cast("ArrayT", weights)


def usable_observations(values: ArrayT, *, name: str = "y") -> ArrayT:
  """Validate a univariate sample and drop the observations that are not one.

  The margins all need the same three steps before they can fit -- require the
  one-dimensional layout, drop ``NaN``, and refuse a column with nothing left
  -- so they share them rather than each spelling the refusal differently.

  Parameters
  ----------
  values : array
      The observations.
  name : str, default="y"
      The argument's name, used in the error messages.

  Returns
  -------
  array
      The finite observations, one dimension, possibly shorter than the input.

  Raises
  ------
  ValueError
      If ``values`` is not one-dimensional, or if no observation survives.
  """
  checked: Any = validate_univariate(values, name=name)
  xp = array_namespace(checked)
  checked = checked[~xp.isnan(checked)]
  if int(checked.shape[0]) == 0:
    raise ValueError(f"{name} has no usable observation")
  return cast("ArrayT", checked)


def reject_covariates(
  part: object, x: Optional[ArrayT], *, name: str = "x"
) -> None:
  """Raise if ``x`` was supplied to a part that fits unconditionally.

  Evaluation ignores covariates a part does not read, since one distribution
  may mix conditional and unconditional parts and they all see the same ``x``.
  Fitting cannot: silently estimating ``f(y)`` when ``f(y | x)`` was asked for
  returns a different model than the caller believes they have.

  Parameters
  ----------
  part : object
      The margin, pair copula or vine being fitted; named in the message.
      Either the instance or the class, so a classmethod may pass ``cls``.
  x : array, shape (n, p), or None, optional
      The covariates the caller passed.
  name : str, default="x"
      The argument's name, used in the error message.

  Returns
  -------
  None

  Raises
  ------
  ValueError
      If ``x`` is not ``None``.
  """
  if x is not None:
    # A classmethod passes `cls`, an instance method `self`; naming
    # `type(cls)` would print the metaclass, which for a Protocol subclass is
    # `_ProtocolMeta` and tells the caller nothing.
    named = part if isinstance(part, type) else type(part)
    raise ValueError(
      f"{named.__name__} does not model covariates, so it cannot be "
      f"fitted with {name}=; write one whose fit reads them and that declares "
      f"supports_covariates, or drop {name}."
    )


def validate_declaration(
  var_type: Optional[str],
  support: Optional[tuple[Optional[float], Optional[float]]],
) -> tuple[Optional[str], Optional[tuple[Optional[float], Optional[float]]]]:
  """Check and normalize what a caller declared about one variable.

  The declaration a margin's ``fit`` / ``select`` / ``from_data`` takes
  keyword-only, checked in one place so every margin refuses the same things
  with the same message.

  Parameters
  ----------
  var_type : {"c", "d", "zi"}, or None, optional
      The variable's type, or ``None`` to leave it to the margin.
  support : tuple of float, or None, optional
      Declared bounds as ``(lo, hi)``, either end ``None`` for unbounded.

  Returns
  -------
  tuple
      The pair, with ``support`` normalized to a 2-tuple.

  Raises
  ------
  ValueError
      If ``var_type`` is not one of the accepted values, or ``support`` is not
      an increasing pair.
  """
  if var_type is not None and var_type not in ("c", "d", "zi"):
    raise ValueError(f"var_type={var_type!r} is not one of ['c', 'd', 'zi']")
  if support is None:
    return var_type, None
  bounds = tuple(support)
  if len(bounds) != 2:
    raise ValueError(f"support must be a (lo, hi) pair; got {bounds!r}")
  lo, hi = bounds
  if lo is not None and hi is not None and not lo < hi:
    raise ValueError(f"support={bounds!r} is not an increasing interval")
  return var_type, (lo, hi)


@contextlib.contextmanager
def extra_required(*, extra: str, requirement: str) -> Iterator[None]:
  """Rewrite an optional dependency's ``ImportError`` to name its extra.

  Guards the ``import`` statement rather than taking a module name to import
  itself. Two reasons: the statement stays where a reader and a type checker
  can see it, and ``importlib.import_module`` would not go through
  ``builtins.__import__`` -- which is how the extras-absent tests simulate an
  absent package, so a guard built on it reports success there.

  Parameters
  ----------
  extra : str
      The extra that installs the dependency, as ``pyvinecopulib[<extra>]``.
  requirement : str
      One sentence saying what needs it, leading the message. It carries its
      own subject and verb, since one class requires a package where a group
      of them require it.

  Yields
  ------
  None
      With the guard installed for the block.

  Raises
  ------
  ImportError
      If the guarded import fails.
  """
  try:
    yield
  except ImportError as e:  # pragma: no cover - exercised in a subprocess
    raise ImportError(
      f"{requirement} Install it with `pip install pyvinecopulib[{extra}]`."
    ) from e
