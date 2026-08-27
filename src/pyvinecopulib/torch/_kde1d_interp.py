"""Pure-torch port of kde1d's cubic interpolation grid.

Reproduces ``kde1d::interp::InterpolationGrid`` (``lib/kde1d/include/kde1d/
interpolation.hpp``) on tensors, so a fitted density can be evaluated on device
and under autograd. Fidelity to the C++ is the contract, including two
deliberate oddities that a "fix" would silently turn into a parity failure:

* **`integrate` has no tail contribution.** It accumulates cell integrals only,
  so the unnormalized integral saturates at the grid's total mass rather than
  approaching one; ``normalize=True`` divides by that mass. The density does
  carry tail mass, so outside the grid the two are not each other's derivative;
  that is upstream's shape, not a slip here.
* **Boundary tangents are zero**, because ``dt0`` / ``dt2`` collapse there,
  which is what makes the extrapolation smooth.

Every function takes ``grid_points`` and ``values`` rather than holding them, so
the caller owns registration and nothing here caches: coefficients and cumulative
cell integrals are recomputed inside the graph on each call. At ``m = 401`` that
is microseconds against an ``n``-sized query, and a cache would both freeze the
gradient and go stale the moment ``values.requires_grad_(True)`` is set --
the opposite call from ``TorchBicop.cache_integrals``, where the cached quantity
is an ``O(m^2)`` two-dimensional integral per query rather than an ``O(m)``
vector shared by the whole batch.
"""

from __future__ import annotations

import torch
from torch import Tensor


def cell_coefs(grid_points: Tensor, values: Tensor) -> Tensor:
  """Monomial coefficients of the cubic on every cell.

  Parameters
  ----------
  grid_points : Tensor, shape (m,)
      Ascending grid.
  values : Tensor, shape (m,)
      Density values on the grid.

  Returns
  -------
  Tensor, shape (m - 1, 4)
      Row ``k`` holds ``(a0, a1, a2, a3)`` for the cell ``[g[k], g[k + 1]]``,
      parameterized on ``[0, 1]``.
  """
  m = grid_points.shape[0]
  k = torch.arange(m - 1, device=grid_points.device)
  k0 = torch.clamp(k - 1, min=0)
  k2 = k + 1
  k3 = torch.clamp(k + 2, max=m - 1)
  dt0 = grid_points[k] - grid_points[k0]
  dt1 = grid_points[k2] - grid_points[k]
  dt2 = grid_points[k3] - grid_points[k2]
  # The tangent is zero where its neighbor collapses. Substituting a dummy
  # denominator before the `where` matters for more than tidiness: dividing by
  # zero and masking the result afterwards still sends a NaN through the
  # backward pass, because the masked branch is multiplied by zero rather than
  # skipped.
  ok0 = dt0 > 0
  ok2 = dt2 > 0
  s0 = torch.where(ok0, dt0, torch.ones_like(dt0))
  s2 = torch.where(ok2, dt2, torch.ones_like(dt2))
  dx1 = torch.where(
    ok0,
    (values[k] - values[k0]) / s0
    - (values[k2] - values[k0]) / (s0 + dt1)
    + (values[k2] - values[k]) / dt1,
    torch.zeros_like(dt0),
  )
  dx2 = torch.where(
    ok2,
    (values[k2] - values[k]) / dt1
    - (values[k3] - values[k]) / (dt1 + s2)
    + (values[k3] - values[k2]) / s2,
    torch.zeros_like(dt2),
  )
  # Rescale to the unit parameterization, then the Schmidt-Hess positivity clip
  # (DOI:10.1007/bf01934097).
  dx1 = torch.maximum(dx1 * dt1, -3.0 * values[k])
  dx2 = torch.minimum(dx2 * dt1, 3.0 * values[k2])
  diff = values[k] - values[k2]
  return torch.stack(
    [
      values[k],
      dx1,
      -3.0 * diff - 2.0 * dx1 - dx2,
      2.0 * diff + dx1 + dx2,
    ],
    dim=-1,
  )


def find_cell(grid_points: Tensor, x: Tensor) -> Tensor:
  """Index of the cell each point falls in.

  The C++ ``find_cell`` bisection: the largest ``k`` with ``g[k] <= x``, clipped
  into ``[0, m - 2]``. The right-hand search side is load-bearing -- at an
  interior knot the cell to the *right* is the one used, which is why the
  batched bicop helper's left-sided search is not interchangeable with this.

  Parameters
  ----------
  grid_points : Tensor, shape (m,)
      Ascending grid.
  x : Tensor, shape (n,)
      Query points.

  Returns
  -------
  Tensor, shape (n,), dtype int64
      Cell indices.
  """
  m = grid_points.shape[0]
  idx = torch.searchsorted(grid_points.contiguous(), x.contiguous(), right=True)
  return torch.clamp(idx - 1, min=0, max=m - 2)


def interpolate(grid_points: Tensor, values: Tensor, x: Tensor) -> Tensor:
  """Evaluate the interpolating spline, with Gaussian-tail extrapolation.

  Parameters
  ----------
  grid_points : Tensor, shape (m,)
      Ascending grid.
  values : Tensor, shape (m,)
      Density values on the grid.
  x : Tensor, shape (n,)
      Query points; ``NaN`` in gives ``NaN`` out.

  Returns
  -------
  Tensor, shape (n,)
      Interpolated values. Not clamped, since the spline can dip below zero
      near a sharp feature and each caller decides what to do about it.
  """
  coefs = cell_coefs(grid_points, values)
  k = find_cell(grid_points, x)
  t = (x - grid_points[k]) / (grid_points[k + 1] - grid_points[k])
  a = coefs[k]
  cubic = _cubic_poly(t, a)
  # Each tail decays from the end it leaves, so both meet the spline
  # continuously: `t` at the left end of the first cell, `t - 1` at the right
  # end of the last.
  out = torch.where(
    t <= 0.0,
    values[k] * torch.exp(-0.5 * t * t),
    torch.where(
      t >= 1.0,
      values[k + 1] * torch.exp(-0.5 * (t - 1.0) * (t - 1.0)),
      cubic,
    ),
  )
  return torch.where(torch.isnan(x), torch.full_like(out, float("nan")), out)


def _cumulative_cell_integrals(grid_points: Tensor, coefs: Tensor) -> Tensor:
  """Mass below each grid point, as a length-``m`` prefix sum."""
  widths = grid_points[1:] - grid_points[:-1]
  full = (
    coefs[:, 0] + coefs[:, 1] / 2.0 + coefs[:, 2] / 3.0 + coefs[:, 3] / 4.0
  ) * widths
  zero = torch.zeros(1, dtype=full.dtype, device=full.device)
  return torch.cat([zero, torch.cumsum(full, dim=0)])


def integrate(
  grid_points: Tensor,
  values: Tensor,
  x: Tensor,
  *,
  normalize: bool = False,
) -> Tensor:
  """Integrate the spline from the grid's left end up to each point.

  The C++ walks the grid once in ascending order of ``x``; the closed form is
  the prefix sum of whole-cell integrals plus the partial cell, which is what
  this evaluates. Below the grid the answer is ``0``; at or above its right end
  it is the total mass, with **no** Gaussian-tail contribution.

  Parameters
  ----------
  grid_points : Tensor, shape (m,)
      Ascending grid.
  values : Tensor, shape (m,)
      Density values on the grid.
  x : Tensor, shape (n,)
      Upper limits; ``NaN`` in gives ``NaN`` out.
  normalize : bool, default=False
      Divide by the total mass, making the result a distribution function.

  Returns
  -------
  Tensor, shape (n,)
      The integral at each point.
  """
  m = grid_points.shape[0]
  coefs = cell_coefs(grid_points, values)
  cum = _cumulative_cell_integrals(grid_points, coefs)
  idx = torch.searchsorted(grid_points.contiguous(), x.contiguous(), right=True)
  # Unlike `find_cell` this admits `m - 1`, the state the C++ loop ends in once
  # the limit is past the grid: every cell has been integrated and there is no
  # partial one left.
  k = torch.clamp(idx - 1, min=0, max=m - 1)
  kc = torch.clamp(k, max=m - 2)
  widths = grid_points[kc + 1] - grid_points[kc]
  t = (x - grid_points[kc]) / widths
  a = coefs[kc]
  partial = (
    a[:, 0] * t
    + a[:, 1] / 2.0 * t**2
    + a[:, 2] / 3.0 * t**3
    + a[:, 3] / 4.0 * t**4
  ) * widths
  res = cum[k] + torch.where(k < m - 1, partial, torch.zeros_like(partial))
  res = torch.where(x <= grid_points[0], torch.zeros_like(res), res)
  res = torch.where(torch.isnan(x), torch.full_like(res, float("nan")), res)
  if normalize:
    res = res / cum[m - 1]
  return res


# The powers are built by repeated multiplication rather than `pow`, and the
# terms summed left to right, because that is what the C++ does: `pow(t, 3)`
# and `t * t * t` need not agree in the last bit, and the inversion below
# iterates on the difference of two of these.
def _cubic_poly(t: Tensor, a: Tensor) -> Tensor:
  """A cell's cubic at position ``t``, one row of coefficients per query."""
  t2 = t * t
  t3 = t2 * t
  return a[:, 0] + a[:, 1] * t + a[:, 2] * t2 + a[:, 3] * t3


def _cubic_indef_integral(t: Tensor, a: Tensor) -> Tensor:
  """Indefinite integral of :func:`_cubic_poly`, zero at ``t = 0``."""
  t2 = t * t
  t3 = t2 * t
  t4 = t3 * t
  return (
    a[:, 0] * t + a[:, 1] / 2.0 * t2 + a[:, 2] / 3.0 * t3 + a[:, 3] / 4.0 * t4
  )


def total_mass(grid_points: Tensor, values: Tensor) -> Tensor:
  """Integral of the spline across the whole grid.

  What :func:`integrate` saturates at, and what ``normalize=True`` divides by.

  Parameters
  ----------
  grid_points : Tensor, shape (m,)
      Ascending grid.
  values : Tensor, shape (m,)
      Density values on the grid.

  Returns
  -------
  Tensor, shape ()
      The total mass.
  """
  return _cumulative_cell_integrals(
    grid_points, cell_coefs(grid_points, values)
  )[-1]


def invert_integral(
  grid_points: Tensor,
  values: Tensor,
  p: Tensor,
  *,
  n_iter: int = 35,
) -> Tensor:
  """Invert the *normalized* integral, cell by cell.

  Reproduces ``kde1d::interp::InterpolationGrid::invert_integral``: locate the
  cell holding the requested cumulative mass, seed a position by interpolating
  that cell's mass linearly, then refine with Newton steps kept inside a
  bracket, bisecting wherever the density is flat. The C++ stops a query as
  soon as its residual is within ``8 eps`` of the total mass; here a query that
  reaches the tolerance is frozen instead, which is the same arithmetic on the
  same iterations. The powers below are built the way the C++ builds them for
  the same reason. What that buys is agreement to a few ULPs rather than an
  equality: the compiled result is not itself portable -- rebuilding kde1d with
  ``-march=native`` and nothing else moves it 19 ULPs -- so no port can be
  bit-identical to every build of it.

  Parameters
  ----------
  grid_points : Tensor, shape (m,)
      Ascending grid.
  values : Tensor, shape (m,)
      Density values on the grid.
  p : Tensor, shape (n,)
      Probabilities. ``p <= 0`` and ``p >= 1`` return the grid's ends exactly.
  n_iter : int, default=35
      Maximum refinement steps per query.

  Returns
  -------
  Tensor, shape (n,)
      The inverse, detached: the iteration carries no usable gradient, and the
      caller reattaches one through the implicit function theorem.
  """
  with torch.no_grad():
    coefs = cell_coefs(grid_points, values)
    cum = _cumulative_cell_integrals(grid_points, coefs)
    total = cum[-1]
    m = grid_points.numel()

    target = p * total
    # `lower_bound` over the interior prefix sums: the first cell whose upper
    # edge already holds the requested mass.
    cell = torch.searchsorted(cum[1:].contiguous(), target.contiguous()).clamp(
      max=m - 2
    )
    lo_g = grid_points[cell]
    width = grid_points[cell + 1] - lo_g
    base = cum[cell]
    cell_mass = cum[cell + 1] - base
    a = coefs[cell]

    wanted = target - base
    empty = ~(cell_mass > 0.0)
    safe_mass = torch.where(empty, torch.ones_like(cell_mass), cell_mass)
    pos = (wanted / safe_mass).clamp(0.0, 1.0)

    lower = torch.zeros_like(pos)
    upper = torch.ones_like(pos)
    tol = 8.0 * torch.finfo(values.dtype).eps * total
    active = ~empty
    for _ in range(n_iter):
      residual = _cubic_indef_integral(pos, a) * width - wanted
      active = active & (residual.abs() > tol)
      if not bool(active.any()):
        break
      lower = torch.where(active & (residual < 0.0), pos, lower)
      upper = torch.where(active & (residual >= 0.0), pos, upper)
      deriv = _cubic_poly(pos, a) * width
      flat = ~(deriv > 0.0)
      step = pos - residual / torch.where(flat, torch.ones_like(deriv), deriv)
      newton = torch.where(flat, lower, step)
      take = ~flat & (newton > lower) & (newton < upper)
      pos = torch.where(
        active, torch.where(take, newton, 0.5 * (lower + upper)), pos
      )

    out = torch.where(empty, grid_points[cell + 1], lo_g + pos * width)
    out = torch.where(p <= 0.0, grid_points[0], out)
    out = torch.where(p >= 1.0, grid_points[-1], out)
    return torch.where(torch.isnan(p), p, out)
