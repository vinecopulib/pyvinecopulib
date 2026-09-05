"""Vine distribution on NumPy.

The concrete :class:`~pyvinecopulib.core.VinedistBase` for the default lane:
NumPy arrays, :class:`~pyvinecopulib.core.Kde1d` margins, and a
:class:`~pyvinecopulib.core.Vinecop`. Everything the class evaluates is
inherited; what it adds is the fitting and the JSON round-trip, both of which
have to name the concrete copula.
"""

from __future__ import annotations

import copy
from typing import Any, ClassVar, Optional, Self

import numpy as np

from ..pyvinecopulib_ext import Kde1d, Vinecop
from .protocols import ControlsLike
from .vinedist_base import VinedistBase

__all__ = ["Vinedist"]


class Vinedist(VinedistBase[np.ndarray]):
  r"""A multivariate distribution built from a vine copula and its margins.

  By Sklar's theorem a joint distribution factorizes into a copula and its
  marginals, and this class is that factorization made into an object:

  .. math::

     F(\mathbf y) = C\bigl(F_1(y_1), \ldots, F_d(y_d)\bigr), \qquad
     \log f(\mathbf y) = \log c\bigl(F_1(y_1), \ldots, F_d(y_d)\bigr)
       + \sum_{j=1}^d \log f_j(y_j).

  The second identity holds for continuous, discrete and mixed margins alike,
  because each margin's ``pdf`` is a density with respect to its own reference
  measure and the copula factor is normalized by the marginal masses. Nothing
  in the evaluation path branches on the variable type.

  This is the NumPy route, so both halves are NumPy ones: a
  :class:`~pyvinecopulib.core.Vinecop` or any other NumPy-side
  :class:`VinecopLike` for the copula, and any :class:`MarginLike` for the
  margins, including distributions from other libraries once passed through
  :func:`pyvinecopulib.margins.as_margin`. A PyTorch part is refused rather
  than silently detached from its graph — use
  :class:`~pyvinecopulib.torch.TorchVinedist` for that.

  :meth:`from_data` fits :class:`~pyvinecopulib.core.Kde1d` margins and a
  :class:`~pyvinecopulib.core.Vinecop`, and :meth:`from_json` reads one back.
  Everything else — the whole evaluation surface — comes from
  :class:`~pyvinecopulib.core.VinedistBase` and is shared with the PyTorch
  route.

  Constructed as ``Vinedist(vinecop, margins)``; the constructor is
  :class:`~pyvinecopulib.core.VinedistBase`'s and is documented there, down to
  what a single broadcast margin means.

  Raises
  ------
  ValueError
      If the number of margins disagrees with the copula's dimension, or if the
      copula's ``var_types`` disagrees with what the margins declare.
  TypeError
      If the copula or a margin comes from :mod:`pyvinecopulib.torch`.

  See Also
  --------
  pyvinecopulib.core.VinedistBase : The array-agnostic base to subclass.
  pyvinecopulib.torch.TorchVinedist : The PyTorch route.
  pyvinecopulib.core.MarginLike : The margin contract.
  pyvinecopulib.core.VinecopLike : The copula contract.
  pyvinecopulib.margins.as_margin : Coerce a foreign distribution to a margin.

  Examples
  --------
  Combine a fitted copula with explicit parametric margins::

      import numpy as np, scipy.stats as st, pyvinecopulib as pv

      u = pv.utils.to_pseudo_obs(y)
      dist = pv.Vinedist(
        pv.Vinecop.from_data(u),
        margins=[st.norm(0, 1), st.gamma(2.0), st.norm(0, 1)],
      )
      dist.logpdf(y)
      dist.sample(100, seeds=[1])
  """

  # The two halves this route fits: NumPy arrays throughout. Plain comments,
  # not `#:` ones -- autosummary cannot resolve an attribute whose value is
  # itself a class, and the declarations are documented on the base.
  vinecop_class: ClassVar[Optional[type]] = Vinecop
  margin_class: ClassVar[Optional[type]] = Kde1d

  #: `_fit_copula` writes weights into a copy of the controls, so both halves
  #: of the fit are weighted by the one argument.
  supports_weighted_copula: bool = True

  def _bind_dist(
    self,
    vinecop: Any,
    margins: Any,
  ) -> None:
    """Install the parts, refusing PyTorch ones.

    The mirror of :class:`~pyvinecopulib.torch.TorchVinedist`'s refusal of a
    NumPy copula, and for the same reason: this class evaluates on NumPy, so a
    torch part would have its gradients detached and would ignore
    ``.to(device)``. The check is by provenance rather than ``isinstance``,
    since :mod:`pyvinecopulib.core` must import without PyTorch.

    Parameters
    ----------
    vinecop : VinecopLike
        The vine copula.
    margins : sequence of MarginLike, or MarginLike
        The margins.

    Returns
    -------
    None

    Raises
    ------
    TypeError
        If any part comes from :mod:`pyvinecopulib.torch`.
    """
    parts = [("vinecop", vinecop)]
    if isinstance(margins, (list, tuple)):
      parts += [(f"margins[{j}]", m) for j, m in enumerate(margins)]
    else:
      parts.append(("margins", margins))
    for name, part in parts:
      module = type(part).__module__
      if module == "torch" or module.startswith(
        ("torch.", "pyvinecopulib.torch")
      ):
        raise TypeError(
          f"Vinedist cannot hold {name}={type(part).__name__}: it evaluates "
          "on NumPy, so a PyTorch part would be detached from its graph and "
          "would ignore `.to(device)`. Use "
          "pyvinecopulib.torch.TorchVinedist instead."
        )
    super()._bind_dist(vinecop, margins)

  @classmethod
  def _coerce_fit_data(
    cls,
    y: Any,
    weights: Optional[Any],
    controls: Optional[ControlsLike],
  ) -> tuple[Any, Any]:
    """Put the fit inputs on NumPy.

    Parameters
    ----------
    y : object
        Observations, including a DataFrame.
    weights : object, or None
        Observation weights.
    controls : ControlsLike, or None
        Unused here; NumPy carries no placement.

    Returns
    -------
    tuple
        The observations as a float array, and the weights unchanged —
        ``validate_weights`` coerces those against the data.
    """
    del controls
    return np.asarray(y, dtype=float), weights

  @classmethod
  def _margin_from_controls(cls, controls: Optional[Any]) -> Any:
    """Build one ``Kde1d`` margin honoring what its controls declare.

    A kernel density takes its variable type and its bounds at construction,
    so a declaration has to reach it before the fit rather than after: an
    unbounded grid is padded past the data, and for a variable that cannot be
    negative that puts mass where nothing can occur.

    Parameters
    ----------
    controls : object, or None
        This variable's marginal configuration.

    Returns
    -------
    Kde1d
        An unfitted margin, bounded and typed as declared.
    """
    from ..margins._resolve import kde_from_controls

    return kde_from_controls(controls)

  @classmethod
  def _fit_copula(
    cls,
    u: Any,
    *,
    var_types: list[str],
    controls: Optional[ControlsLike],
    structure: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> Any:
    """Fit a ``Vinecop`` on the pseudo-observations.

    Parameters
    ----------
    u : array, shape (n, d + k), dtype float
        The copula-scale layout the margins produced.
    var_types : list of str
        One ``"c"`` or ``"d"`` per variable.
    controls : FitControlsVinecop, or None
        Copula fit controls.
    structure : RVineStructure, or None, optional
        A fixed structure; selected from the data when ``None``.
    weights : array, shape (n,), or None, optional
        Observation weights, written into a copy of the controls so both halves
        of the fit are weighted by the one argument.

    Returns
    -------
    Vinecop
        The fitted copula.
    """
    from ..pyvinecopulib_ext import FitControlsVinecop

    resolved: Any = FitControlsVinecop() if controls is None else controls
    if weights is not None:
      # The explicit argument governs both halves. Copy first so overriding a
      # controls object's weights cannot mutate an object the caller still owns.
      resolved = copy.deepcopy(resolved)
      resolved.weights = np.asarray(weights, dtype=float)
    return super()._fit_copula(
      u,
      var_types=var_types,
      controls=resolved,
      structure=structure,
      weights=None,
    )

  @classmethod
  def from_json(cls, json: str) -> Self:
    """Instantiate from a JSON string.

    Parameters
    ----------
    json : str
        A string produced by :meth:`to_json`.

    Returns
    -------
    Vinedist
        The deserialized distribution, with a ``Vinecop`` copula.

    Raises
    ------
    ValueError
        If the payload's version is unrecognized, or a margin's ``kind`` has no
        registered reader.
    """
    from ..pyvinecopulib_ext import Vinecop
    from ._serialization import (
      MARGIN_JSON_VERSION,
      loads,
      margin_from_json,
    )

    payload = loads(json)
    if payload.get("version") != MARGIN_JSON_VERSION:
      raise ValueError(
        f"unsupported Vinedist JSON version {payload.get('version')!r}; this "
        f"build reads version {MARGIN_JSON_VERSION}"
      )
    return cls(
      Vinecop.from_json(payload["copula"]),
      [margin_from_json(m) for m in payload["margins"]],
    )
