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

from array_api_compat import array_namespace

from .protocols import ArrayT

__all__ = [
  "reject_covariates",
  "usable_observations",
  "validate_weights",
]


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
  # Every entry being a drop marker leaves nothing to fit: `Kde1d` rescales by
  # the surviving sum, so it reaches the caller as a crash rather than an
  # error, and `utils.wdm` refuses the same input for the same reason.
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


def reject_array_controls(part: object, controls: object) -> None:
  """Raise if an array landed in the ``controls`` slot.

  Every estimator in the package takes the observations, then ``controls``.
  ``Kde1d`` is the documented exception -- its second positional argument is
  ``weights`` -- so ``kde.fit(x, w)`` is a spelling a reader carries over, and
  on any other margin it binds the weights to ``controls``, where they are
  ignored: an unweighted fit under a weighted-looking call.

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
    "weights are the keyword-only `weights=`; `Kde1d` is the one class whose "
    "second positional argument is `weights`."
  )


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
