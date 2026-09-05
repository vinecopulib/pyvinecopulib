"""Fit configuration for a univariate margin."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Sequence

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

  The marginal counterpart of ``FitControlsBicop``, and the second half of the
  specification a vine distribution takes: ``margins=`` says *which* margin
  class each variable gets, and this says *how* to fit or select it. Both are
  resolved per variable by the same rule, so one call may bound the two
  variables whose bounds are known and leave the rest alone
  (:func:`~pyvinecopulib.margins.resolve_margin_controls`).

  Two fields configure a family search and are read by
  :meth:`~pyvinecopulib.margins.SciPyMargin.select`; the other two carry
  what the *caller* knows about the variable and the margin cannot infer. A
  margin whose type or support is fixed at construction keeps what it was
  built with -- the declaration is a default, not an instruction -- but the
  library honors it when it is the one constructing the margin, which is what
  makes a bounded default reachable without naming a class.

  Attributes
  ----------
  family_set : sequence of str, or None, default=None
      Candidate families for a search, named as in the ecosystem the margin
      belongs to (``"gamma"``, ``"lognorm"``). ``None`` searches the curated
      set admissible for the variable, which is the recommendation: a blind
      sweep ranks a family that misstates its own support above the truth.
  selection_criterion : {"aic", "bic", "aicc"}, default="aic"
      ``"aic"`` (the default), ``"bic"`` or ``"aicc"``.
  var_type : {"c", "d", "zi"}, or None, default=None
      ``"c"``, ``"d"`` or ``"zi"``. Worth setting whenever the caller knows
      the variable type and the sample does not show it -- a count column
      whose smallest observation is 3 reads as continuous.
  support : tuple of float, or None, default=None
      Declared bounds as ``(lo, hi)``, either of which may be ``None`` for
      unbounded on that side. What a bound *means* differs by variable type;
      see the ``concepts-kde-margins`` section of the concepts page.
  on_failure : {"raise", "fallback"}, default="raise"
      ``"raise"`` (the default) reports every candidate and why it lost;
      ``"fallback"`` substitutes a kernel-density margin with one warning.

  Raises
  ------
  ValueError
      If ``selection_criterion``, ``on_failure`` or ``var_type`` is not one of
      the accepted values, if ``family_set`` is empty, or if ``support`` is
      not an ordered pair.

  See Also
  --------
  pyvinecopulib.margins.resolve_margin_controls : Expand one per variable.

  Examples
  --------
  >>> from pyvinecopulib.margins import FitControlsMargin
  >>> FitControlsMargin(support=(0.0, None)).support
  (0.0, None)
  """

  family_set: Optional[Sequence[str]] = None
  selection_criterion: str = "aic"
  var_type: Optional[str] = None
  support: Optional[tuple[Optional[float], Optional[float]]] = None
  on_failure: str = "raise"

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
    if self.var_type is not None and self.var_type not in ("c", "d", "zi"):
      raise ValueError(
        f"var_type={self.var_type!r} is not one of ['c', 'd', 'zi']"
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
    if self.support is not None:
      bounds = tuple(self.support)
      if len(bounds) != 2:
        raise ValueError(f"support must be a (lo, hi) pair; got {bounds!r}")
      lo, hi = bounds
      if lo is not None and hi is not None and not lo < hi:
        raise ValueError(f"support={bounds!r} is not an increasing interval")
      self.support = (lo, hi)

  def to_dict(self) -> dict[str, Any]:
    """Settings as a dictionary, which is what makes this a ``ControlsLike``.

    Returns
    -------
    dict
        One entry per setting.
    """
    return {
      "family_set": None if self.family_set is None else list(self.family_set),
      "selection_criterion": self.selection_criterion,
      "var_type": self.var_type,
      "support": self.support,
      "on_failure": self.on_failure,
    }
