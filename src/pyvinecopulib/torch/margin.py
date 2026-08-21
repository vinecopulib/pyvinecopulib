"""PyTorch univariate margin built on a ``torch.distributions`` family.

A ``torch.distributions.Distribution`` is not a :class:`torch.nn.Module`: it has
no ``.to(device)``, and held as a plain attribute it contributes nothing to a
module's ``state_dict``. :class:`TorchMargin` therefore registers the
*parameters* — as ``torch.nn.Parameter`` or as buffers — and rebuilds the
distribution on every call, which is the same shape
:class:`~pyvinecopulib.torch.TorchBicop` uses for its interpolation grid.

See Also
--------
pyvinecopulib.core.MarginBase : The contract this fills in.
TorchVinedist : The joint distribution these margins go into.
"""

from __future__ import annotations

from itertools import chain
from typing import Any, Callable, Iterable, Mapping, Optional, Union, cast

import torch
from torch import Tensor
from torch.distributions import Distribution

from ..core import MarginBase

__all__ = ["TorchMargin"]

#: Accepted spellings of the ``parameters`` argument: anything ``dict`` takes.
ParameterSpec = Union[Mapping[str, Any], Iterable[tuple[str, Any]]]


def _bound_to_float(bound: Any, unbounded: float) -> float:
  """Read one support endpoint as a plain float.

  A parameter-derived endpoint (``Uniform``'s, say) arrives as a tensor that may
  carry a gradient, and converting one to a Python scalar without detaching it
  first is a warning waiting to happen.

  Parameters
  ----------
  bound : object
      The endpoint as the constraint spells it: a float, a tensor, or ``None``.
  unbounded : float
      What to return when the constraint does not pin this side.

  Returns
  -------
  float
      The endpoint.
  """
  if bound is None:
    return unbounded
  return float(bound.detach() if isinstance(bound, Tensor) else bound)


def _support_bounds(distribution: Distribution) -> tuple[float, float]:
  """Read ``(lo, hi)`` off a distribution's support constraint.

  Parameters
  ----------
  distribution : torch.distributions.Distribution
      The distribution to inspect.

  Returns
  -------
  tuple of float
      The bounds, unbounded on either side the constraint does not pin.
  """
  constraint = distribution.support
  return (
    _bound_to_float(getattr(constraint, "lower_bound", None), float("-inf")),
    _bound_to_float(getattr(constraint, "upper_bound", None), float("inf")),
  )


def _implements(distribution: Distribution, name: str) -> bool:
  """Whether a family overrides one of ``Distribution``'s stubs.

  ``Distribution.cdf`` and ``Distribution.icdf`` exist on every distribution and
  raise :exc:`NotImplementedError` unless the concrete family overrides them, so
  member names carry no information here. Comparing the resolved function
  against the base one does, and costs no evaluation.

  Parameters
  ----------
  distribution : torch.distributions.Distribution
      The distribution to inspect.
  name : str
      Method name, ``"cdf"`` or ``"icdf"``.

  Returns
  -------
  bool
      ``True`` if the family supplies its own implementation.
  """
  return getattr(type(distribution), name, None) is not getattr(
    Distribution, name
  )


class TorchMargin(MarginBase[Tensor], torch.nn.Module):
  """Univariate margin on a ``torch.distributions`` family.

  A :class:`~pyvinecopulib.core.MarginBase` that is also a
  :class:`torch.nn.Module`, so its parameters move with ``.to(device)``, appear
  in ``state_dict()``, and are visible to an optimizer. The distribution object
  itself is rebuilt from those parameters on every call and never stored.

  Parameters
  ----------
  factory : callable
      Called with the registered parameters as keywords to produce the
      distribution. Normally the family itself, e.g.
      ``torch.distributions.Normal``.
  parameters : mapping or iterable of pair, optional
      Parameter values, keyed by the keyword ``factory`` expects. Anything
      :class:`dict` accepts; the default is no parameters at all.
  trainable : bool, default=True
      Register floating-point parameters as ``torch.nn.Parameter`` rather
      than as buffers. Integer parameters are always buffers, since a tensor of
      integers cannot carry a gradient.
  validate_args : bool or None, default=None
      Forwarded to ``factory``; ``None`` leaves the family's own default, which
      is what makes an optimizer step that leaves the parameter domain raise
      rather than return silent ``nan``.
  device : torch.device or None, default=None
      Placement of the registered tensors.
  dtype : torch.dtype, default=torch.float64
      Precision of the registered floating-point tensors, matching the rest of
      :mod:`pyvinecopulib.torch`.

  Raises
  ------
  NotImplementedError
      If the family's support is discrete, or if the family does not implement
      ``cdf``.

  See Also
  --------
  from_distribution : Lift an already-constructed distribution.
  pyvinecopulib.core.MarginBase : The contract this fills in.
  pyvinecopulib.margins.as_margin : NumPy-side coercion of foreign objects.

  Notes
  -----
  **Coverage.** ``cdf`` is what puts a variable on the copula scale, so a family
  without one cannot be a margin at all and is rejected at construction. Of the
  families in ``torch.distributions`` (checked against torch 2.13):

  - ``cdf`` and ``icdf`` both: ``Normal``, ``Uniform``, ``Exponential``,
    ``Cauchy``, ``Laplace``, ``LogNormal``, ``Weibull``, ``Gumbel``,
    ``HalfNormal``.
  - ``cdf`` only: ``Gamma``, ``Chi2``. :meth:`icdf` inverts ``cdf`` by
    bisection over :attr:`support`, so these work too.
  - neither: ``Beta``, ``StudentT``, ``FisherSnedecor``, ``VonMises``. No
    fallback can rescue them — inverting a ``cdf`` that does not exist is not a
    matter of effort — so they raise. Reach for
    :class:`pyvinecopulib.core.Kde1d`, or a SciPy distribution through
    :class:`pyvinecopulib.core.Vinedist` on NumPy.
  - discrete (``Poisson``, ``Binomial``, ``Categorical``, …): rejected. The
    torch lane is continuous-only for now; the copula needs a left-limit ``cdf``
    at every atom, which the discrete cascade on
    :class:`pyvinecopulib.core.Vinedist` provides and this one does not.

  **Autograd.** ``Normal.cdf`` is differentiable with respect to both its
  parameters and its argument. ``Gamma.cdf`` is differentiable with respect to
  its argument only — ``igamma`` has no derivative for its first argument — so
  a Gamma margin can be fitted by its own marginal likelihood (``log_prob`` *is*
  differentiable in the concentration) but cannot be optimized end to end
  through a copula, which needs ``dF/dtheta``. That asymmetry is a reason to
  keep two-step estimation the default even here.

  **Constrained parameters.** Values are stored on the family's natural scale,
  unconstrained. An optimizer is free to step a scale negative; wrap the
  parameter with ``torch.distributions.transform_to`` if that matters, or
  pass ``trainable=False`` to freeze it.

  Examples
  --------
  A standard-normal margin whose location is learnable::

      import torch
      from pyvinecopulib.torch import TorchMargin

      margin = TorchMargin(
        torch.distributions.Normal,
        {"loc": 0.0, "scale": 1.0},
      )
      margin.cdf(torch.tensor([-1.0, 0.0, 1.0], dtype=torch.float64))
      list(margin.state_dict())  # ['loc', 'scale']
  """

  def __init__(
    self,
    factory: Callable[..., Distribution],
    parameters: ParameterSpec = (),
    *,
    trainable: bool = True,
    validate_args: Optional[bool] = None,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> None:
    # Initialize nn.Module explicitly: TorchMargin also subclasses MarginBase
    # (a Protocol-derived ABC), whose __init__ chain would otherwise shadow
    # nn.Module's under super().
    torch.nn.Module.__init__(self)
    values = dict(parameters)
    self._factory = factory
    self._validate_args = validate_args
    self._parameter_names = tuple(values)
    self._fallback_dtype = dtype
    self._fallback_device = device
    self.family_name = getattr(factory, "__name__", type(factory).__name__)

    for name, value in values.items():
      # Probe the kind, then convert from the source at the target precision in
      # one step: a Python float goes through torch's default float32 first
      # otherwise, and widening a rounded float32 back to float64 leaves the
      # parameter wrong in its eighth digit. Integers keep their own dtype.
      probe = torch.as_tensor(value)
      target = dtype if probe.is_floating_point() else probe.dtype
      tensor = (
        torch.as_tensor(value, dtype=target, device=device).detach().clone()
      )
      if trainable and tensor.is_floating_point():
        self.register_parameter(name, torch.nn.Parameter(tensor))
      else:
        self.register_buffer(name, tensor)

    probe = self.distribution
    if getattr(probe.support, "is_discrete", False):
      raise NotImplementedError(
        f"{self.family_name} is a discrete family, and pyvinecopulib.torch is "
        "continuous-only: the copula needs a left-limit cdf at every atom. Use "
        "pyvinecopulib.core.Vinedist with pyvinecopulib.margins for discrete "
        "or mixed data."
      )
    if not _implements(probe, "cdf"):
      raise NotImplementedError(
        f"{self.family_name} does not implement cdf, so it cannot serve as a "
        "margin: the copula scale is F(x). Nothing can be inverted or "
        "differenced around that. Use pyvinecopulib.core.Kde1d, or a "
        "SciPy distribution through pyvinecopulib.core.Vinedist."
      )
    self._has_icdf = _implements(probe, "icdf")

  # --------------------------------------------------------------------- #
  # Constructor                                                            #
  # --------------------------------------------------------------------- #

  @classmethod
  def from_distribution(
    cls,
    distribution: Distribution,
    *,
    trainable: bool = True,
    device: Optional[torch.device] = None,
    dtype: torch.dtype = torch.float64,
  ) -> "TorchMargin":
    """Lift an already-constructed ``torch.distributions`` object.

    The family and its parameter names are read off the object — the latter
    from ``arg_constraints``, which is how ``torch.distributions`` itself
    declares them — and the values are copied into registered tensors.

    Parameters
    ----------
    distribution : torch.distributions.Distribution
        The distribution to lift. Its parameters are copied, not shared, so
        mutating it afterwards does not affect the margin.
    trainable : bool, default=True
        See the class docstring.
    device : torch.device or None, default=None
        Placement of the registered tensors; ``None`` keeps the originals'.
    dtype : torch.dtype, default=torch.float64
        Precision of the registered floating-point tensors.

    Returns
    -------
    TorchMargin
        A margin equivalent to ``distribution``.

    Raises
    ------
    ValueError
        If a name in ``arg_constraints`` is not readable off the object.

    Notes
    -----
    A family parameterized by ``probs`` *or* ``logits`` declares both and
    refuses to be handed both; the one it was built with is the one lifted.
    """
    names = list(distribution.arg_constraints)
    if "probs" in names and "logits" in names:
      # One parameter spelled two ways: the family derives whichever it was not
      # given as a `lazy_property`, and refuses to be handed both. Keep the one
      # already materialized on the object, which is the one it was built with.
      names.remove("probs" if "logits" in vars(distribution) else "logits")
    missing = [name for name in names if not hasattr(distribution, name)]
    if missing:
      raise ValueError(
        f"{type(distribution).__name__} declares {missing} in "
        "arg_constraints but does not expose them; pass the parameters to "
        "TorchMargin(...) explicitly"
      )
    return cls(
      type(distribution),
      {name: getattr(distribution, name) for name in names},
      trainable=trainable,
      device=device,
      dtype=dtype,
    )

  # --------------------------------------------------------------------- #
  # Structure                                                              #
  # --------------------------------------------------------------------- #

  @property
  def distribution(self) -> Distribution:
    """The distribution, rebuilt from the currently registered parameters.

    A fresh object every call, so an in-place parameter update, an optimizer
    step, or a ``.to(device)`` is picked up immediately.

    Returns
    -------
    torch.distributions.Distribution
        The distribution ``factory`` produces for the current parameters.
    """
    kwargs: dict[str, Any] = {
      name: getattr(self, name) for name in self._parameter_names
    }
    if self._validate_args is not None:
      kwargs["validate_args"] = self._validate_args
    return self._factory(**kwargs)

  @property
  def parameter_names(self) -> tuple[str, ...]:
    """Names of the registered parameters, in the order given.

    Spelled ``parameter_names`` rather than ``parameters``, which
    :class:`torch.nn.Module` already uses for the tensors themselves.

    Returns
    -------
    tuple of str
        The keywords ``factory`` is called with.
    """
    return self._parameter_names

  @property
  def support(self) -> tuple[float, float]:
    """Closed bounds of the support.

    Read off the distribution rather than cached, since a family such as
    ``Uniform`` derives its support from its parameters.

    Returns
    -------
    tuple of float
        ``(lo, hi)``, infinite where unbounded.
    """
    return _support_bounds(self.distribution)

  # --------------------------------------------------------------------- #
  # Evaluation                                                             #
  # --------------------------------------------------------------------- #

  def pdf(self, y: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Density of the margin.

    Parameters
    ----------
    y : Tensor, shape (n,), dtype float
        Observations on the original scale.
    x : Tensor, shape (n, k), or None, optional
        Ignored; this margin is unconditional.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Density values.
    """
    return self.distribution.log_prob(y).exp()

  def logpdf(self, y: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Log-density of the margin.

    Taken straight from the family's ``log_prob`` rather than through
    :meth:`pdf`, which keeps the far tails finite.

    Parameters
    ----------
    y : Tensor, shape (n,), dtype float
        Observations on the original scale.
    x : Tensor, shape (n, k), or None, optional
        Ignored; this margin is unconditional.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Log-density values.
    """
    return self.distribution.log_prob(y)

  def log_prob(self, y: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Alias of :meth:`logpdf`, the ``torch.distributions`` spelling.

    Parameters
    ----------
    y : Tensor, shape (n,), dtype float
        Observations on the original scale.
    x : Tensor, shape (n, k), or None, optional
        Ignored; this margin is unconditional.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Log-density values.
    """
    return self.logpdf(y, x=x)

  def cdf(self, y: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Distribution function of the margin.

    Parameters
    ----------
    y : Tensor, shape (n,), dtype float
        Observations on the original scale.
    x : Tensor, shape (n, k), or None, optional
        Ignored; this margin is unconditional.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Distribution values in ``[0, 1]``.
    """
    return self.distribution.cdf(y)

  def icdf(self, p: Tensor, *, x: Optional[Tensor] = None) -> Tensor:
    """Inverse distribution function of the margin.

    Uses the family's own ``icdf`` where it has one, and otherwise inverts
    :meth:`cdf` by bisection over :attr:`support` — which is what makes
    ``Gamma`` and ``Chi2`` usable.

    The bisection runs outside autograd, as
    :meth:`pyvinecopulib.torch.TorchVinecop.sample` does: the derivative of a
    fixed number of bisection steps is not the derivative of the inverse, and
    the families that need the fallback have no ``dF/dtheta`` to propagate
    anyway. The result of the closed-form branch stays differentiable.

    Parameters
    ----------
    p : Tensor, shape (n,), dtype float
        Probabilities in ``[0, 1]``.
    x : Tensor, shape (n, k), or None, optional
        Ignored; this margin is unconditional.

    Returns
    -------
    Tensor, shape (n,), dtype float
        Quantiles on the original scale.
    """
    if self._has_icdf:
      return self.distribution.icdf(p)
    with torch.no_grad():
      return super().icdf(p)

  def _ref_tensor(self) -> Tensor:
    """A registered tensor to crib dtype/device from."""
    for tensor in chain(self.parameters(), self.buffers()):
      return tensor
    # A factory closing over its own parameters registers none, and so is not
    # placed anywhere; fall back to what construction asked for.
    return torch.empty(
      0, dtype=self._fallback_dtype, device=self._fallback_device
    )

  def _sample_uniform(self, n: int, seeds: list[int]) -> Tensor:
    """Draw ``n`` uniforms on the registered tensors' dtype/device."""
    ref = self._ref_tensor()
    generator: Optional[torch.Generator] = None
    if seeds:
      generator = torch.Generator(device=ref.device).manual_seed(int(seeds[0]))
    return torch.rand(
      n, generator=generator, dtype=ref.dtype, device=ref.device
    )

  def __repr__(self) -> str:
    """Return a structural representation of the margin.

    Returns
    -------
    str
        The family and its current parameter values.
    """
    body = ", ".join(
      f"{name}={cast(Tensor, getattr(self, name)).detach().cpu().tolist()}"
      for name in self._parameter_names
    )
    return f"{type(self).__name__}({self.family_name}({body}))"
