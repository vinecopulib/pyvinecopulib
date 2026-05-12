"""Model confidence set (MCS) survivor selection for vine forests.

Private helpers used by :class:`VineForestBase` to prune the random-search
candidate pool down to a model confidence set after fitting. Two
selectors are exposed:

- :func:`da_mcs_marg` — single-split DA test with marginal
  per-model error control (Hansen, Lunde & Nason 2011; Kim, Olsen,
  Nagler & Vatter 2025).
- :func:`da_mcs_unif` — dual-split variant with uniform (familywise)
  error control via an adaptive pre-screening step.

Both consume an ``(n, p)`` matrix of per-sample losses across ``p``
candidate models and return a boolean survivor mask plus the
underlying test statistics. See the supplement to Kim et al. (2025)
for the theory.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm


def argmin_except_self(mu: np.ndarray) -> np.ndarray:
  """Return, for each i, the argmin of ``mu`` excluding index i.

  Used to pair each candidate model with its best-performing
  competitor. Exploits the fact that excluding self only matters when
  i actually attains the minimum: every other row gets the global min;
  the global-min row gets the second min.
  """
  p = mu.size
  v1 = mu.min()
  mins = np.flatnonzero(mu == v1)  # sorted ascending
  k = np.empty(p, dtype=int)

  if mins.size >= 2:
    # Two (or more) rows attain the minimum: each can use any other
    # minimum-attaining row as its partner.
    j1, j2 = mins[0], mins[1]
    k.fill(j1)
    k[j1] = j2
    return k

  # Unique minimum: mask it and find the runner-up.
  j1 = mins[0]
  tmp = mu.copy()
  tmp[j1] = np.inf
  j2 = int(tmp.argmin())

  k.fill(j1)
  k[j1] = j2
  return k


def pairwise_sd_for_pairs(
  X: np.ndarray, mu: np.ndarray, k: np.ndarray, ddof: int = 1
) -> np.ndarray:
  """SD of column differences X[:, i] - X[:, k[i]] across all i."""
  n, p = X.shape
  r = n - ddof
  C = X - mu
  var = (C * C).sum(axis=0) / r

  # `k` from argmin_except_self has at most two distinct values.
  partners = np.unique(k)
  assert partners.size in (1, 2), "argmin_except_self() guarantees at most two"

  var_diff = np.empty(p)
  for partner in partners:
    idx = np.flatnonzero(k == partner)
    cov_partner = (C * C[:, [partner]]).sum(axis=0) / r
    var_diff[idx] = var[idx] + var[partner] - 2.0 * cov_partner[idx]

  var_diff = np.maximum(var_diff, 0.0)
  return np.sqrt(var_diff)


def pairwise_sd_for_pairs_einsum(
  X: np.ndarray, mu: np.ndarray, k: np.ndarray, ddof: int = 1
) -> np.ndarray:
  """einsum-based variant of :func:`pairwise_sd_for_pairs`.

  Avoids the (n, p)-sized centring buffer; faster for tall ``X``.
  """
  n, p = X.shape
  r = n - ddof

  ss = np.einsum("ij,ij->j", X, X)
  var = (ss - n * mu * mu) / r

  partners = np.unique(k)
  assert 1 <= partners.size <= 2

  var_diff = np.empty(p)

  cov_by_partner: dict[int, np.ndarray] = {}
  for partner in partners:
    cov_by_partner[int(partner)] = (
      np.einsum("ij,i->j", X, X[:, partner]) - n * mu * mu[partner]
    ) / r

  for partner in partners:
    idx = np.flatnonzero(k == partner)
    cov = cov_by_partner[int(partner)][idx]
    var_diff[idx] = var[idx] + var[partner] - 2.0 * cov

  var_diff = np.maximum(var_diff, 0.0)
  return np.sqrt(var_diff)


def da_mcs_unif(
  data: np.ndarray,
  alpha: float = 0.05,
  use_einsum: bool = True,
  randomize: bool = True,
  rng: np.random.Generator | None = None,
) -> dict[str, np.ndarray]:
  r"""Dual-split DA test with uniform error control.

  Splits the sample into halves ``(ind1, ind2)``, then splits ``ind1``
  again into ``(ind11, ind12)``. The first inner split pre-screens to
  estimate the survivor-set size :math:`\widetilde M`, and the outer
  split runs the final DA test with a Bonferroni-corrected threshold
  :math:`\Phi^{-1}(1 - \alpha / \widetilde M)`. This controls the
  familywise probability of excluding any optimal model uniformly.

  Parameters
  ----------
  data : ndarray of shape (n, p)
      Per-sample loss for each of the ``p`` candidate models.
  alpha : float, default=0.05
      Significance level.
  use_einsum : bool, default=True
      Use the einsum-based standard-deviation helper.
  randomize : bool, default=True
      Sample the splits uniformly. If ``False``, splits are
      deterministic ``[0, m), [m, n)`` halves.
  rng : numpy.random.Generator, optional
      Random generator; only used when ``randomize=True``.

  Returns
  -------
  dict
      ``{"decision": bool ndarray (p,), "stats": float ndarray (p,),
      "stats_pre": float ndarray (p,)}``.
  """
  n, p = data.shape
  m = n // 2
  m1 = m // 2

  if randomize:
    if rng is None:
      rng = np.random.default_rng()
    ind1 = rng.choice(n, size=m, replace=False)
    ind2 = np.setdiff1d(np.arange(n), ind1)
    ind11 = rng.choice(ind1, size=m1, replace=False)
    ind12 = np.setdiff1d(ind1, ind11)
  else:
    ind1 = np.arange(m)
    ind2 = np.arange(m, n)
    ind11 = np.arange(m1)
    ind12 = np.arange(m1, m)

  X1, X2 = data[ind1, :], data[ind2, :]
  X11, X12 = data[ind11, :], data[ind12, :]

  mu1 = X1.mean(axis=0)
  mu2 = X2.mean(axis=0)
  mu11 = X11.mean(axis=0)
  mu12 = X12.mean(axis=0)

  # Pre-screen on (X11, X12) at the relaxed level n_val^{-1/2}.
  k_pre = argmin_except_self(mu11)
  mu_pre = mu12 - mu12[k_pre]
  sd_fn = pairwise_sd_for_pairs_einsum if use_einsum else pairwise_sd_for_pairs
  sd_pre = sd_fn(X12, mu12, k_pre)
  stat_pre = np.sqrt(m1) * (mu_pre / sd_pre)
  thr_pre = norm.ppf(1 - 1 / np.sqrt(m))
  decision_pre = stat_pre < thr_pre
  S = int(decision_pre.sum()) or 1

  # Final test on (X1, X2) with Bonferroni correction.
  k_fin = argmin_except_self(mu1)
  mu_fin = mu2 - mu2[k_fin]
  sd_fin = sd_fn(X2, mu2, k_fin)
  stat_fin = np.sqrt(m) * (mu_fin / sd_fin)
  thr_fin = norm.ppf(1 - alpha / S)
  decision = stat_fin < thr_fin

  return {"decision": decision, "stats": stat_fin, "stats_pre": stat_pre}


def da_mcs_marg(
  data: np.ndarray,
  alpha: float = 0.05,
  use_einsum: bool = True,
  randomize: bool = True,
  rng: np.random.Generator | None = None,
) -> dict[str, np.ndarray]:
  r"""Single-split DA test with marginal per-model error control.

  Splits the sample into halves once, then computes the
  pairwise-difference test statistic at the un-corrected level
  :math:`\Phi^{-1}(1 - \alpha)`. The resulting survivor mask has
  marginal (per-model) coverage but no familywise correction;
  cheaper than :func:`da_mcs_unif` and often the default choice
  for forest-style ensembles.

  Parameters
  ----------
  data : ndarray of shape (n, p)
      Per-sample loss for each of the ``p`` candidate models.
  alpha : float, default=0.05
      Per-model significance level.
  use_einsum : bool, default=True
      Use the einsum-based standard-deviation helper.
  randomize : bool, default=True
      Sample the split uniformly. If ``False``, the split is the
      deterministic halves ``[0, m), [m, n)``.
  rng : numpy.random.Generator, optional
      Random generator; only used when ``randomize=True``.

  Returns
  -------
  dict
      ``{"decision": bool ndarray (p,), "stats": float ndarray (p,),
      "threshold": float scalar}``.
  """
  n, p = data.shape
  m = n // 2

  if randomize:
    if rng is None:
      rng = np.random.default_rng()
    ind1 = rng.choice(n, size=m, replace=False)
    ind2 = np.setdiff1d(np.arange(n), ind1)
  else:
    ind1 = np.arange(m)
    ind2 = np.arange(m, n)

  X1, X2 = data[ind1, :], data[ind2, :]

  mu1 = X1.mean(axis=0)
  mu2 = X2.mean(axis=0)

  k = argmin_except_self(mu1)
  mu = mu2 - mu2[k]
  sd_fn = pairwise_sd_for_pairs_einsum if use_einsum else pairwise_sd_for_pairs
  sd = sd_fn(X2, mu2, k)
  stat = np.sqrt(m) * (mu / sd)
  thr = np.asarray(norm.ppf(1 - alpha))
  decision = stat < thr

  return {"decision": decision, "stats": stat, "threshold": thr}
