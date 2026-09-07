"""The package's validators, so one rule and one message serve every layer.

Every check a margin, pair copula, vine or vine distribution performs on its
inputs lives here: the univariate and covariate layouts, observation weights,
and the fit-time refusal of covariates a part cannot read. Each takes ``name=``
so the message names the argument the caller actually passed.
"""

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
  if weights.ndim != 1 or int(weights.shape[0]) != int(values.shape[0]):
    raise ValueError(
      f"{name} must have shape ({int(values.shape[0])},), with one weight "
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
  # Every entry being a drop marker leaves nothing to fit: `Kde1d` rescales by
  # the surviving sum, so it reaches the caller as a crash rather than an
  # error, and `utils.wdm` refuses the same input for the same reason.
  kept = weights[~xp.isnan(weights)]
  if not float(xp.sum(kept)) > 0.0:
    raise ValueError(
      f"{name} must leave at least one observation standing; NaN and 0 mark a "
      "dropped observation, and every entry is one"
    )
  return weights


def usable_observations(values: Any, *, name: str = "y") -> Any:
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
  values = validate_univariate(values, name=name)
  xp = array_namespace(values)
  values = values[~xp.isnan(values)]
  if int(values.shape[0]) == 0:
    raise ValueError(f"{name} has no usable observation")
  return values


def reject_covariates(part: Any, x: Optional[Any], *, name: str = "x") -> None:
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
  x : array, shape (n, p), or None
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


def reject_array_controls(part: Any, controls: Any) -> None:
  """Raise if an array landed in the ``controls`` slot.

  Every estimator in the package takes the observations, then ``controls``.
  The compiled ``Kde1d`` is the documented exception -- its second positional
  argument is ``weights`` -- so ``kde.fit(x, w)`` is a spelling a reader
  carries over, and on any other margin it binds the weights to ``controls``,
  where they are ignored: an unweighted fit under a weighted-looking call.

  Nothing in the library passes an array here, so refusing one costs nothing
  and turns that typo into a message naming the keyword to use.

  Parameters
  ----------
  part : object
      The margin being fitted; named in the message. Either the instance or
      the class, so a classmethod may pass ``cls``.
  controls : object
      Whatever arrived in the controls slot.

  Returns
  -------
  None

  Raises
  ------
  TypeError
      If ``controls`` looks like an array rather than a configuration object.
  """
  if controls is None or hasattr(controls, "to_dict"):
    return
  if not any(hasattr(controls, name) for name in ("shape", "__array__")):
    return
  named = part if isinstance(part, type) else type(part)
  raise TypeError(
    f"{named.__name__} received an array where `controls` goes. Observation "
    "weights are the keyword-only `weights=`; the compiled `Kde1d` is the one "
    "class whose second positional argument is `weights`."
  )
