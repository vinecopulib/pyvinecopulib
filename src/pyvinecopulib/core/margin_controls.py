"""Fit configuration for a univariate margin."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, fields
from typing import Any

__all__ = ["FitControlsMargin"]

#: Information criteria a family search may be scored on. All three are
#: minimized, and all three charge for freely estimated parameters.
CRITERIA: tuple[str, ...] = ("aic", "bic", "aicc")

#: What a family search does when no candidate is admissible. ``"raise"``
#: reports every family and its cause; ``"fallback"`` substitutes a
#: kernel-density margin and warns once.
ON_FAILURE: tuple[str, ...] = ("raise", "fallback")


@dataclass
class FitControlsMargin:
  """Configuration for fitting or selecting one margin.

  The marginal counterpart of ``FitControlsBicop``, and the marginal half of a
  vine distribution's fit: the distribution's ``margin_class`` says *which*
  margin each variable gets, and this says *how* to fit or select it. It is
  resolved per variable four ways -- one object broadcast, a length-``d``
  sequence, or a mapping keyed by position or name -- so one call may configure
  the variables that need it and leave the rest alone
  (:func:`~pyvinecopulib.margins.resolve_margin_controls`).

  Fit configuration only: every field here is about *how* to estimate the
  margin. What the caller knows about the variable -- its type and its bounds
  -- is a *declaration*, and travels as the keyword-only ``var_type`` and
  ``support`` arguments of ``fit`` / ``select`` / ``from_data``, exactly as
  ``var_types`` does on ``Bicop.from_data``.

  ``weights`` is here for the same reason it is on ``FitControlsBicop``:
  observation weights are a setting, so they travel in the controls at every
  level. There is no rule propagating them from one controls object to
  another, which is what lets a vine distribution's margins and its copula be
  weighted differently -- or one of them weighted and the other not.

  Attributes
  ----------
  family_set : sequence of str, or None, default=None
      Candidate families for a search, named as in the ecosystem the margin
      belongs to (``"gamma"``, ``"lognorm"``). ``None`` searches the curated
      set admissible for the variable, which is the recommendation: an
      unfiltered sweep ranks a family that misstates its own support above
      the truth.
  selection_criterion : {"aic", "bic", "aicc"}, default="aic"
      ``"aic"`` (the default), ``"bic"`` or ``"aicc"``.
  on_failure : {"raise", "fallback"}, default="raise"
      ``"raise"`` (the default) reports every candidate and why it lost;
      ``"fallback"`` substitutes a kernel-density margin with one warning.
  weights : array, shape (n,), or None, optional
      Observation weights. Refused by a margin that declares
      ``supports_weights`` ``False``, rather than fitted away.

  Raises
  ------
  ValueError
      If ``selection_criterion`` or ``on_failure`` is not one of the accepted
      values, if ``family_set`` is empty, or if ``weights`` is empty.

  See Also
  --------
  pyvinecopulib.margins.resolve_margin_controls : Expand one per variable.

  Examples
  --------
  >>> from pyvinecopulib.margins import FitControlsMargin
  >>> FitControlsMargin(family_set=["gamma", "lognorm"]).family_set
  ['gamma', 'lognorm']
  """

  family_set: Sequence[str] | None = None
  selection_criterion: str = "aic"
  on_failure: str = "raise"
  # `Any` rather than `Array`: the weights are handed on to the
  # `ArrayT`-parameterized `validate_weights`, which a bare `Array` does not
  # satisfy, and this class is not itself generic -- one controls type serves
  # a NumPy margin and a torch one alike.
  weights: Any | None = None

  def __post_init__(self) -> None:
    """Validate the settings.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If any setting is outside its accepted set.
    """
    if self.selection_criterion not in CRITERIA:
      raise ValueError(
        f"selection_criterion={self.selection_criterion!r} is not one of "
        f"{list(CRITERIA)}"
      )
    if self.on_failure not in ON_FAILURE:
      raise ValueError(
        f"on_failure={self.on_failure!r} is not one of {list(ON_FAILURE)}"
      )
    if self.family_set is not None:
      families = list(self.family_set)
      if not families:
        raise ValueError(
          "family_set is empty, so no candidate could be selected; pass None "
          "for the curated set"
        )
      if not all(isinstance(f, str) for f in families):
        raise ValueError("family_set must name families as strings")
      self.family_set = families
    if self.weights is not None and len(self.weights) == 0:
      raise ValueError("weights is empty; pass None for an unweighted fit")

  def to_dict(self) -> dict[str, Any]:
    """Settings as a dictionary, which is what makes this a ``ControlsLike``.

    Returns
    -------
    dict
        One entry per setting, keyed by the field name.
    """
    return {f.name: getattr(self, f.name) for f in fields(self)}
