"""PyTorch univariate margins on a ``torch.distributions`` family.

``TorchMargin`` is the marginal half a ``TorchVinedist`` composes with its
copula: a ``MarginBase`` that is also a ``torch.nn.Module``, so a family's
parameters move with ``.to(device)``, appear in ``state_dict()`` and are
visible to an optimizer.

A ``torch.distributions.Distribution`` is not a module — it has no
``.to(device)``, and held as a plain attribute it contributes nothing to a
``state_dict`` — so what a margin registers is the *parameters*, and the
distribution is rebuilt from them on every call. That is the same shape
``TorchBicop`` uses for its interpolation grid.

Continuous families only: a margin with atoms needs a left-limit ``cdf``,
which ``torch.distributions`` does not expose. ``TorchKde1d`` is the discrete
and zero-inflated margin on this lane; ``Vinedist`` with
``pyvinecopulib.margins`` is the NumPy one.

See Also
--------
pyvinecopulib.core.MarginBase : The contract this fills in.
TorchKde1d : The margin for discrete and zero-inflated variables.
TorchVinedist : The joint distribution these margins go into.
"""

from __future__ import annotations

from itertools import chain
from typing import Any, Callable, Iterable, Mapping, Optional, Union, cast

import torch
from torch import Tensor
from torch.distributions import Distribution

from ..core import MarginBase
from ..core.margin_base import support_of

__all__ = ["TorchMargin"]

#: Accepted spellings of the ``parameters`` argument: anything ``dict`` takes.
ParameterSpec = Union[Mapping[str, Any], Iterable[tuple[str, Any]]]


def _implements(distribution: Distribution, name: str) -> bool:
  """Whether a family supplies its own ``cdf`` or ``icdf``.

  ``Distribution`` defines both and raises ``NotImplementedError`` unless the
  concrete family overrides them, so the member is there either way and its
  presence says nothing; which function the name resolves to does.

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

  A ``MarginBase`` that is also a ``torch.nn.Module``: the parameter values are
  registered tensors, so they move with ``.to(device)``, appear in
  ``state_dict()`` and are visible to an optimizer, while the distribution
  itself is rebuilt from them on every call and never stored.

  What a caller supplies is the family and its parameters — as ``factory`` plus
  ``parameters`` here, or as an already-built distribution through
  ``from_distribution``. Everything else reads off those:

  - ``distribution`` is the family at the parameter values as they stand now,
    which is what every member below evaluates.
  - ``parameter_names`` names the parameters in the order they were given, and
    ``support`` reports the bounds those values imply.
  - ``pdf`` / ``logpdf`` (``log_prob`` is the ``torch.distributions`` spelling
    of the same method) / ``cdf`` / ``icdf`` are the evaluation surface, and
    ``MarginBase`` derives ``cdf_left`` / ``loglik`` / ``sample`` from them.

  The parameters are the whole specification, so a margin is ready to evaluate
  as soon as it is constructed and the inherited ``fit`` / ``select`` /
  ``from_data`` estimators have nothing to estimate and raise: hand the
  constructor an estimate, or let an optimizer move the registered parameters.
  ``TorchKde1d`` is the margin on this lane that fits itself to data.

  Parameters
  ----------
  factory : callable
      Called with the registered parameters as keywords to produce the
      distribution. Normally the family itself, e.g.
      ``torch.distributions.Normal``.
  parameters : mapping or iterable of pair, optional
      Parameter values, keyed by the keyword ``factory`` expects. Anything
      ``dict`` accepts; the default is no parameters at all, for a factory that
      closes over its own.
  trainable : bool, default=True
      Register floating-point parameters as ``torch.nn.Parameter`` rather
      than as buffers. Integer parameters are always buffers, since a tensor of
      integers cannot carry a gradient.
  validate_args : bool or None, default=None
      Forwarded to ``factory``; ``None`` leaves the family's own default, which
      normally checks its parameters, so a value outside the family's domain
      raises rather than returning silent ``nan``.
  device : torch.device or None, default=None
      Placement of the registered tensors.
  dtype : torch.dtype, default=torch.float64
      Precision of the registered floating-point tensors, matching the rest of
      ``pyvinecopulib.torch``.

  Raises
  ------
  NotImplementedError
      If the family is discrete, or if it does not implement ``cdf``.

  See Also
  --------
  from_distribution : Lift an already-constructed distribution.
  pyvinecopulib.core.MarginBase : The contract this fills in.
  pyvinecopulib.margins.as_margin : NumPy-side coercion of foreign objects.

  Notes
  -----
  **Coverage.** ``cdf`` is what puts a variable on the copula scale, so a family
  without one cannot be a margin at all and is rejected at construction.
  Grouping the families of ``torch.distributions`` by what they implement
  (torch 2.13):

  - ``cdf`` and ``icdf`` both: ``Normal``, ``Uniform``, ``Exponential``,
    ``Cauchy``, ``HalfCauchy``, ``Laplace``, ``LogNormal``, ``Weibull``,
    ``Gumbel``, ``HalfNormal``, ``Pareto``, ``InverseGamma``, ``Kumaraswamy``.
  - ``cdf`` only: ``Gamma``, ``Chi2``. ``icdf`` inverts ``cdf`` by
    bisection over ``support``, so these work too.
  - neither: ``Beta``, ``StudentT``, ``FisherSnedecor``, ``VonMises``. No
    fallback can rescue them — inverting a ``cdf`` that does not exist is not a
    matter of effort — so they raise. Reach for ``Kde1d``, or a SciPy
    distribution through ``Vinedist`` on NumPy.
  - discrete (``Poisson``, ``Binomial``, ``Categorical``, …): rejected. A
    margin with atoms needs a left-limit ``cdf`` at each of them, which none of
    them exposes; ``TorchKde1d``, or ``Vinedist`` with
    ``pyvinecopulib.margins``, covers discrete and mixed data instead.

  **Autograd.** ``Normal.cdf`` is differentiable with respect to its argument
  and to both parameters. ``Gamma.cdf`` is differentiable with respect to its
  argument and to ``rate``, but not to ``concentration`` — ``igamma`` has no
  derivative for its first argument — so a Gamma margin's shape can be
  estimated from its own marginal likelihood (``log_prob`` *is* differentiable
  in it) but not end to end through a copula, which reaches a margin as
  ``F(y)`` alone. That asymmetry is a reason to keep two-step estimation the
  default even here.

  **Constrained parameters.** Values are stored on the family's natural scale,
  unconstrained. An optimizer is free to step a scale negative; wrap the
  parameter with ``torch.distributions.transform_to`` if that matters, or
  pass ``trainable=False`` to freeze it.

  Examples
  --------
  A standard-normal margin whose parameters are learnable::

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

    The family and the names of its parameters are read off the object, and the
    values are copied into registered tensors.

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
        If the family declares a parameter in ``arg_constraints`` that the
        object does not expose; pass the parameters to the constructor
        instead.

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
    """The family at the currently registered parameter values.

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
    ``torch.nn.Module`` already uses for the tensors themselves.

    Returns
    -------
    tuple of str
        The keywords ``factory`` is called with.
    """
    return self._parameter_names

  @property
  def support(self) -> tuple[float, float]:
    """Closed bounds of the support.

    Follows the current parameter values, so a family such as ``Uniform``,
    whose support is set by them, stays right after an update. ``icdf``
    brackets its search here.

    Returns
    -------
    tuple of float
        ``(lo, hi)``, infinite where unbounded.
    """
    return support_of(self.distribution)

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

    Finite far into the tails, where ``pdf`` itself underflows to zero.

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
    """The ``torch.distributions`` spelling of ``logpdf``.

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

    The family's own ``icdf`` where it has one, and otherwise ``cdf`` inverted
    numerically over ``support`` — which is what makes ``Gamma`` and ``Chi2``
    usable as margins. Only the first of those carries a gradient: the
    derivative of a fixed number of bisection steps is not the derivative of an
    inverse, so the numerical branch returns a detached quantile, and with it a
    detached ``sample``.

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

  def _apply(self, fn: Any, *args: Any, **kwargs: Any) -> "TorchMargin":
    """Keep the fallback placement current across a ``.to()``.

    A factory that closes over its own parameters registers none, so the
    fallback dtype and device are the only record of where this margin lives,
    and they are where ``sample`` draws.
    """
    out = super()._apply(fn, *args, **kwargs)
    probe = fn(
      torch.empty(0, dtype=self._fallback_dtype, device=self._fallback_device)
    )
    self._fallback_dtype = probe.dtype
    self._fallback_device = probe.device
    return out

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
