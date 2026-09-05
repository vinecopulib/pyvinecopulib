"""Shared validation for row-aligned conditional and weighted evaluations."""

from __future__ import annotations

from typing import Any, Optional

from array_api_compat import array_namespace

__all__: list[str] = []


def validate_univariate(values: Any, *, name: str = "y") -> Any:
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
  x: Optional[Any], n_rows: int, *, name: str = "x"
) -> None:
  """Require a two-dimensional covariate matrix row-aligned with the data.

  Parameters
  ----------
  x : array, or None
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
    raise ValueError(
      f"{name} must have shape (n, p), with one row per observation; "
      f"got {shape}"
    )
  rows = int(x.shape[0])
  if rows != n_rows:
    raise ValueError(
      f"{name} must have one row per observation; got {rows} rows for "
      f"{n_rows} observations"
    )


def validate_weights(
  weights: Optional[Any], values: Any, *, name: str = "weights"
) -> Optional[Any]:
  """Normalize and validate one real, finite, nonnegative weight per row.

  Parameters
  ----------
  weights : array, or None
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

  Raises
  ------
  ValueError
      If the weights are not one per observation, not finite, or negative.
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
  if weights.ndim != 1 or int(weights.shape[0]) != int(values.shape[0]):
    raise ValueError(
      f"{name} must have shape ({int(values.shape[0])},), with one weight "
      f"per observation; got {tuple(weights.shape)}"
    )
  if not xp.isdtype(weights.dtype, ("real floating", "integral")):
    raise TypeError(
      f"{name} must have a real numeric dtype; got {weights.dtype}"
    )
  if not bool(xp.all(xp.isfinite(weights))):
    raise ValueError(f"{name} must contain only finite values")
  if bool(xp.any(weights < 0)):
    raise ValueError(f"{name} must be nonnegative")
  return weights
