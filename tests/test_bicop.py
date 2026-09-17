import os
from collections.abc import Callable

import numpy as np
import pytest

import pyvinecopulib as pv


@pytest.fixture
def gumbel() -> pv.Bicop:
  """A one-parameter Gumbel, the fixture most of these assertions read."""
  cop = pv.Bicop(family=pv.families.gumbel, rotation=90)
  cop.rotation = 0
  cop.parameters = np.array([[3.0]])
  return cop


#: Two rows, enough for every per-observation member to answer with a shape.
_U = np.array([[0.1, 0.2], [0.3, 0.4]])


def test_a_default_bicop_is_the_independence_copula() -> None:
  """And so is one fitted to two rows that carry no dependence."""
  for cop in (pv.Bicop(), pv.Bicop.from_data(_U, pv.FitControlsBicop())):
    assert cop.family == pv.families.indep
    assert cop.rotation == 0
    assert cop.parameters.shape == (0, 0)
    assert cop.var_types == ["c", "c"]


def test_it_round_trips_through_json_and_both_file_formats(
  unique_json_path: str,
) -> None:
  """Text for a `.json` name, binary for `.cbor`, and the model survives both."""
  cop = pv.Bicop.from_data(_U, pv.FitControlsBicop())

  def same(other: pv.Bicop) -> None:
    assert cop.family == other.family
    assert cop.rotation == other.rotation
    assert cop.parameters.shape == other.parameters.shape
    assert cop.var_types == other.var_types

  same(pv.Bicop.from_json(cop.to_json()))

  filename = os.fspath(unique_json_path)
  cop.to_file(filename)
  assert os.path.exists(filename)
  same(pv.Bicop.from_file(filename))
  # A non-`.cbor` name keeps writing JSON text.
  with open(filename, "rb") as f:
    assert f.read(1) == b"{"

  cbor_filename = filename.replace(".json", ".cbor")
  cop.to_file(cbor_filename)
  same(pv.Bicop.from_file(cbor_filename))
  with open(cbor_filename, "rb") as f:
    assert f.read(1) != b"{"


def test_the_writable_properties_refuse_what_they_cannot_hold(
  gumbel: pv.Bicop,
) -> None:
  """`rotation` takes the four right angles; `parameters` and `var_types` set."""
  assert gumbel.rotation == 0
  with pytest.raises(RuntimeError):
    gumbel.rotation = 45

  assert gumbel.parameters.shape == (1, 1)
  assert gumbel.parameters[0, 0] == 3.0

  gumbel.var_types = ["d", "d"]
  assert gumbel.var_types == ["d", "d"]


def test_the_read_only_properties_are_read_only(gumbel: pv.Bicop) -> None:
  """`npars` is the fit's own count and no caller's to set."""
  assert isinstance(gumbel.tau, float)
  assert gumbel.npars == 1
  with pytest.raises(AttributeError):
    # Through `setattr`, so `ty` does not reject the assignment whose
    # runtime refusal is the claim.
    setattr(gumbel, "npars", 2)  # noqa: B010


def test_every_evaluation_member_answers_with_a_shape(
  gumbel: pv.Bicop,
) -> None:
  """Including a single row, which used to come back a scalar (#169, #170)."""
  one = gumbel.pdf(np.array([[0.1, 0.2]]))
  assert isinstance(one, np.ndarray) and one.shape == (1,)

  for method in ("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"):
    values = getattr(gumbel, method)(_U)
    assert isinstance(values, np.ndarray)
    assert values.shape == (2,)


def test_the_fit_criteria_are_floats(gumbel: pv.Bicop) -> None:
  """`loglik` and the three criteria the compiled class has always had."""
  assert isinstance(gumbel.loglik(_U), float)
  assert isinstance(gumbel.aic(_U), float)
  assert isinstance(gumbel.bic(_U), float)
  assert isinstance(gumbel.mbic(_U, 0.9), float)


def test_it_renders(gumbel: pv.Bicop) -> None:
  """`repr` and `str` both answer, which the notebooks depend on."""
  assert isinstance(repr(gumbel), str)
  assert isinstance(str(gumbel), str)


def test_the_tau_maps_require_the_family_s_own_parameter_shape(
  gumbel: pv.Bicop,
) -> None:
  """A wrong shape used to read past the family's own storage and answer."""
  assert isinstance(gumbel.parameters_to_tau(gumbel.parameters), float)
  with pytest.raises(RuntimeError, match="wrong shape"):
    gumbel.parameters_to_tau(np.array([[0.5, 0.6], [0.7, 0.8]]))
  assert isinstance(gumbel.tau_to_parameters(0.5), np.ndarray)


def test_the_parameter_bounds_are_the_family_s(gumbel: pv.Bicop) -> None:
  """Gumbel's parameter runs from 1 to 50."""
  assert gumbel.parameters_lower_bounds == np.array([1.0])
  assert gumbel.parameters_upper_bounds == np.array([50.0])


def test_sampling_and_the_two_fitting_verbs(gumbel: pv.Bicop) -> None:
  """`sample` answers `(n, 2)`, and `fit` / `select` both run."""
  assert gumbel.sample(100, False, seeds=[]).shape == (100, 2)
  gumbel.fit(_U, pv.FitControlsBicop())
  gumbel.select(_U, pv.FitControlsBicop())


_PER_ROW_METHODS = ["pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"]


def _per_row_params(
  cop: pv.Bicop, n: int, rng: np.random.RandomState
) -> np.ndarray:
  """Build an (n, p) matrix of valid per-row parameters around ``cop``'s."""
  base = cop.parameters.ravel()
  p = base.shape[0]
  return base[None, :] * (1.0 + 0.1 * rng.uniform(-1.0, 1.0, size=(n, p)))


@pytest.mark.parametrize(
  ("family", "parameters", "rotation"),
  [
    (pv.families.clayton, np.array([[2.0]]), 0),
    (pv.families.clayton, np.array([[2.5]]), 90),
    (pv.families.bb1, np.array([[1.0], [1.5]]), 0),
  ],
)
def test_bicop_per_row_parameters(
  family: pv.BicopFamily, parameters: np.ndarray, rotation: int
) -> None:
  rng = np.random.RandomState(42)
  n = 50
  u = rng.uniform(0.05, 0.95, size=(n, 2))
  cop = pv.Bicop(family=family, rotation=rotation, parameters=parameters)
  p = cop.parameters.shape[0]
  pars = _per_row_params(cop, n, rng)
  pars_const = np.tile(cop.parameters.ravel(), (n, 1))

  for method in _PER_ROW_METHODS:
    fn = getattr(cop, method)
    # Shape: (n, 2) data + (n, p) parameters -> (n,) output.
    values = fn(u, parameters=pars)
    assert isinstance(values, np.ndarray)
    assert values.shape == (n,)

    # Identity parity: constant per-row parameters equal to the object's own
    # parameters reproduce the state-based (single-argument) result.
    np.testing.assert_allclose(
      fn(u, parameters=pars_const), fn(u), rtol=1e-9, atol=1e-12
    )

    # Threading must not change the result.
    np.testing.assert_allclose(
      fn(u, parameters=pars, num_threads=2), values, rtol=1e-12, atol=1e-14
    )

    # Strong per-row parity: row i with parameters pars[i] equals a fresh
    # copula built with those parameters evaluated at row i.
    for i in (0, n // 2, n - 1):
      single = pv.Bicop(
        family=family, rotation=rotation, parameters=pars[i].reshape(p, 1)
      )
      np.testing.assert_allclose(
        values[i], getattr(single, method)(u[i : i + 1])[0], rtol=1e-9
      )

  # loglik with per-row parameters: scalar, NaN-ignoring sum of log-densities,
  # and matches the state-based loglik under constant parameters.
  ll = cop.loglik(u, parameters=pars)
  assert isinstance(ll, float)
  np.testing.assert_allclose(
    ll, np.nansum(np.log(cop.pdf(u, parameters=pars))), rtol=1e-9
  )
  np.testing.assert_allclose(
    cop.loglik(u, parameters=pars_const), cop.loglik(u), rtol=1e-9, atol=1e-12
  )


def test_bicop_per_row_parameters_errors() -> None:
  rng = np.random.RandomState(0)
  n = 20
  u = rng.uniform(0.05, 0.95, size=(n, 2))
  cop = pv.Bicop(family=pv.families.clayton, parameters=np.array([[2.0]]))
  pars = _per_row_params(cop, n, rng)

  # Nonparametric families are not supported.
  data = rng.uniform(size=(200, 2))
  tll = pv.Bicop.from_data(
    data, pv.FitControlsBicop(family_set=[pv.families.tll])
  )
  with pytest.raises(RuntimeError):
    tll.pdf(u, parameters=np.ones((n, 1)))

  # Wrong number of parameter rows (must equal u.rows()).
  with pytest.raises(RuntimeError):
    cop.pdf(u, parameters=pars[:-1])

  # Wrong number of parameter columns (must equal the family's parameter count).
  with pytest.raises(RuntimeError):
    cop.pdf(u, parameters=np.ones((n, 2)))

  # Non-finite parameters are rejected.
  pars_nan = pars.copy()
  pars_nan[0, 0] = np.nan
  with pytest.raises(RuntimeError):
    cop.pdf(u, parameters=pars_nan)

  # Out-of-bounds parameters are rejected.
  pars_oob = pars.copy()
  pars_oob[0, 0] = -5.0
  with pytest.raises(RuntimeError):
    cop.pdf(u, parameters=pars_oob)


_FIRST_ORDER_DERIVS = [
  "pdf_deriv",
  "hfunc1_deriv",
  "hfunc2_deriv",
  "logpdf_deriv",
]
_SECOND_ORDER_DERIVS = [
  "pdf_deriv2",
  "hfunc1_deriv2",
  "hfunc2_deriv2",
  "logpdf_deriv2",
]


@pytest.mark.parametrize(
  ("family", "parameters", "rotation"),
  [
    (pv.families.clayton, np.array([[2.0]]), 0),
    (pv.families.clayton, np.array([[2.5]]), 90),
    (pv.families.gumbel, np.array([[2.0]]), 0),
    (pv.families.gaussian, np.array([[0.5]]), 0),
    (pv.families.student, np.array([[0.5], [4.0]]), 0),
    (pv.families.bb1, np.array([[1.0], [1.5]]), 0),
  ],
)
def test_bicop_deriv_matches_finite_differences(
  family: pv.BicopFamily, parameters: np.ndarray, rotation: int
) -> None:
  rng = np.random.RandomState(42)
  n = 50
  u = rng.uniform(0.1, 0.9, size=(n, 2))
  cop = pv.Bicop(family=family, rotation=rotation, parameters=parameters)
  p = cop.parameters.shape[0]

  # First derivative w.r.t. each parameter vs central finite differences of
  # fresh copulas at theta +/- h (derivatives are w.r.t. the natural, positive
  # parameters of the rotated copula, so the same rotation is kept).
  for k in range(p):
    h = 1e-5 * max(1.0, abs(cop.parameters[k, 0]))
    pars_up, pars_lo = cop.parameters.copy(), cop.parameters.copy()
    pars_up[k, 0] += h
    pars_lo[k, 0] -= h
    up = pv.Bicop(family=family, rotation=rotation, parameters=pars_up)
    lo = pv.Bicop(family=family, rotation=rotation, parameters=pars_lo)
    fd = (up.pdf(u) - lo.pdf(u)) / (2 * h)
    np.testing.assert_allclose(
      cop.pdf_deriv(u, f"par{k + 1}"), fd, rtol=1e-4, atol=1e-6
    )

  # First derivative w.r.t. the arguments vs central finite differences in u.
  for j, sel in enumerate(("u1", "u2")):
    h = 1e-6
    u_up, u_lo = u.copy(), u.copy()
    u_up[:, j] += h
    u_lo[:, j] -= h
    fd = (cop.pdf(u_up) - cop.pdf(u_lo)) / (2 * h)
    np.testing.assert_allclose(cop.pdf_deriv(u, sel), fd, rtol=1e-4, atol=1e-5)

  # Quotient identity: d log c = (d c) / c.
  np.testing.assert_allclose(
    cop.logpdf_deriv(u, "par1"),
    cop.pdf_deriv(u, "par1") / cop.pdf(u),
    rtol=1e-8,
    atol=1e-10,
  )

  # h-function shortcuts: differentiating an h-function w.r.t. its
  # non-conditioning argument yields the copula density.
  np.testing.assert_allclose(
    cop.hfunc1_deriv(u, "u2"), cop.pdf(u), rtol=1e-10, atol=1e-12
  )
  np.testing.assert_allclose(
    cop.hfunc2_deriv(u, "u1"), cop.pdf(u), rtol=1e-10, atol=1e-12
  )

  # Second derivative selectors are order-invariant, and a mixed selector
  # matches finite differences of the analytic first derivative.
  np.testing.assert_allclose(
    cop.pdf_deriv2(u, "par1u1"),
    cop.pdf_deriv2(u, "u1par1"),
    rtol=1e-12,
    atol=1e-14,
  )
  h = 1e-6
  u_up, u_lo = u.copy(), u.copy()
  u_up[:, 0] += h
  u_lo[:, 0] -= h
  fd = (cop.pdf_deriv(u_up, "par1") - cop.pdf_deriv(u_lo, "par1")) / (2 * h)
  np.testing.assert_allclose(
    cop.pdf_deriv2(u, "par1u1"), fd, rtol=1e-3, atol=1e-4
  )

  # Shorthands: "par" means "par1"; a single component in a second-order
  # selector means differentiating twice.
  np.testing.assert_allclose(
    cop.pdf_deriv(u, "par"), cop.pdf_deriv(u, "par1"), rtol=1e-12
  )
  np.testing.assert_allclose(
    cop.pdf_deriv2(u, "par"), cop.pdf_deriv2(u, "par1par1"), rtol=1e-12
  )
  np.testing.assert_allclose(
    cop.pdf_deriv2(u, "u1"), cop.pdf_deriv2(u, "u1u1"), rtol=1e-12
  )


def test_bicop_deriv_per_row_parameters() -> None:
  rng = np.random.RandomState(42)
  n = 50
  u = rng.uniform(0.1, 0.9, size=(n, 2))
  cop = pv.Bicop(family=pv.families.clayton, parameters=np.array([[2.0]]))
  pars = _per_row_params(cop, n, rng)
  pars_const = np.tile(cop.parameters.ravel(), (n, 1))
  selectors = dict.fromkeys(_FIRST_ORDER_DERIVS, "par1") | dict.fromkeys(
    _SECOND_ORDER_DERIVS, "par1u1"
  )

  for method, sel in selectors.items():
    fn = getattr(cop, method)
    values = fn(u, sel, parameters=pars)
    assert isinstance(values, np.ndarray) and values.shape == (n,)

    # Constant per-row parameters reproduce the state-based call.
    np.testing.assert_allclose(
      fn(u, sel, parameters=pars_const), fn(u, sel), rtol=1e-9, atol=1e-12
    )

    # Threading must not change the result.
    np.testing.assert_allclose(
      fn(u, sel, parameters=pars, num_threads=2), values, rtol=1e-12, atol=1e-14
    )

    # Strong per-row parity vs fresh copulas built with the row's parameters.
    for i in (0, n // 2, n - 1):
      single = pv.Bicop(
        family=pv.families.clayton, parameters=pars[i].reshape(1, 1)
      )
      np.testing.assert_allclose(
        values[i], getattr(single, method)(u[i : i + 1], sel)[0], rtol=1e-9
      )


def test_bicop_deriv_errors() -> None:
  rng = np.random.RandomState(0)
  u = rng.uniform(0.1, 0.9, size=(20, 2))
  cop = pv.Bicop(family=pv.families.clayton, parameters=np.array([[2.0]]))

  # Empty, malformed, and out-of-range selectors.
  for bad in ("", "foo", "par2", "u3"):
    with pytest.raises(RuntimeError):
      cop.pdf_deriv(u, bad)

  # A first-order method rejects a two-component selector.
  with pytest.raises(RuntimeError):
    cop.pdf_deriv(u, "par1u1")

  # A second-order method rejects a three-component selector.
  with pytest.raises(RuntimeError):
    cop.pdf_deriv2(u, "par1u1u2")

  # Derivatives require continuous variable types.
  discrete = pv.Bicop(
    family=pv.families.clayton,
    parameters=np.array([[2.0]]),
    var_types=["d", "d"],
  )
  with pytest.raises(RuntimeError):
    discrete.pdf_deriv(u, "par1")


@pytest.mark.xfail(
  reason=(
    "The pinned lib/vinecopulib branch (feature/rvine-trees-perf) drops the "
    "KernelBicop analytic argument-derivative override (vinecopulib#694) and "
    "rejects derivatives for nonparametric copulas, so Bicop.pdf_deriv raises "
    "for tll. Pending an upstream decision on whether tll argument gradients "
    "are restored; un-xfail if they are."
  ),
  strict=False,
  raises=RuntimeError,
)
def test_bicop_deriv_tll() -> None:
  rng = np.random.RandomState(42)
  data = pv.to_pseudo_obs(rng.normal(size=(500, 2)))
  tll = pv.Bicop.from_data(
    data, controls=pv.FitControlsBicop(family_set=[pv.families.tll])
  )
  u = rng.uniform(0.2, 0.8, size=(50, 2))

  # Argument gradients work: exact slope of the bilinear interpolation grid,
  # cross-checked against central finite differences (the surface is piecewise
  # linear in each argument, so away from grid knots FD is near-exact).
  for j, sel in enumerate(("u1", "u2")):
    grad = tll.pdf_deriv(u, sel)
    assert np.isfinite(grad).all()
    h = 1e-7
    u_up, u_lo = u.copy(), u.copy()
    u_up[:, j] += h
    u_lo[:, j] -= h
    fd = (tll.pdf(u_up) - tll.pdf(u_lo)) / (2 * h)
    np.testing.assert_allclose(grad, fd, rtol=1e-3, atol=1e-3)

  # Quotient identity for the log-density gradient.
  np.testing.assert_allclose(
    tll.logpdf_deriv(u, "u1"),
    tll.pdf_deriv(u, "u1") / tll.pdf(u),
    rtol=1e-8,
    atol=1e-10,
  )

  # The non-conditioning-argument shortcut still works (equals the density).
  np.testing.assert_allclose(
    tll.hfunc1_deriv(u, "u2"), tll.pdf(u), rtol=1e-10, atol=1e-12
  )

  # Parameter selectors, second-order derivatives, and actual h-function
  # derivatives are undefined for tll.
  with pytest.raises(RuntimeError):
    tll.pdf_deriv(u, "par1")
  with pytest.raises(RuntimeError):
    tll.pdf_deriv2(u, "u1")
  with pytest.raises(RuntimeError):
    tll.hfunc1_deriv(u, "u1")


def test_families_analytic_derivs() -> None:
  group = pv.families.analytic_derivs
  assert isinstance(group, list)
  assert all(isinstance(f, pv.BicopFamily) for f in group)
  # Currently every parametric family has closed forms (vinecopulib#683/#687);
  # assert bb1/tawn membership explicitly so an upstream change is caught
  # for that reason.
  assert set(group) == set(pv.families.parametric)
  assert pv.families.bb1 in group
  assert pv.families.tawn in group
  assert pv.families.tll not in group


def test_bicop_taildep_and_beta() -> None:
  # Clayton theta=2: lower tail dependence 2^(-1/theta), no upper.
  clayton = pv.Bicop(family=pv.families.clayton, parameters=np.array([[2.0]]))
  td = clayton.taildep
  assert isinstance(td, np.ndarray) and td.shape == (2, 2)
  np.testing.assert_allclose(td[0, 0], 2 ** (-1 / 2), rtol=1e-12)
  assert td[1, 1] == 0.0 and td[0, 1] == 0.0 and td[1, 0] == 0.0

  # Gumbel theta=2: upper tail dependence 2 - 2^(1/theta), no lower.
  gumbel = pv.Bicop(family=pv.families.gumbel, parameters=np.array([[2.0]]))
  np.testing.assert_allclose(gumbel.taildep[1, 1], 2 - 2**0.5, rtol=1e-12)
  assert gumbel.taildep[0, 0] == 0.0

  # Gaussian: no tail dependence in any corner.
  gaussian = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
  np.testing.assert_allclose(gaussian.taildep, np.zeros((2, 2)), atol=1e-12)

  # Student t: all four corners positive; the concordant (diagonal) corners
  # are equal, and so are the discordant ones (rho replaced by -rho).
  student = pv.Bicop(
    family=pv.families.student, parameters=np.array([[0.5], [4.0]])
  )
  td = student.taildep
  assert (td > 0).all()
  np.testing.assert_allclose(td[0, 0], td[1, 1], rtol=1e-12)
  np.testing.assert_allclose(td[0, 1], td[1, 0], rtol=1e-12)
  assert td[0, 0] > td[0, 1]

  # Rotations: 180 degrees swaps lower/upper; 90 degrees moves the dependence
  # to the off-diagonal corners and flips beta's sign.
  clayton180 = pv.Bicop(
    family=pv.families.clayton, rotation=180, parameters=np.array([[2.0]])
  )
  np.testing.assert_allclose(
    clayton180.taildep[1, 1], clayton.taildep[0, 0], rtol=1e-12
  )
  assert clayton180.taildep[0, 0] == 0.0
  clayton90 = pv.Bicop(
    family=pv.families.clayton, rotation=90, parameters=np.array([[2.0]])
  )
  td90 = clayton90.taildep
  assert td90[0, 0] == 0.0 and td90[1, 1] == 0.0
  assert td90[0, 1] + td90[1, 0] > 0
  assert clayton90.beta < 0 < clayton.beta

  # Blomqvist's beta identity: beta = 4 * C(0.5, 0.5) - 1.
  for cop in (clayton, gumbel, gaussian, student, clayton90, clayton180):
    np.testing.assert_allclose(
      cop.beta,
      4 * cop.cdf(np.array([[0.5, 0.5]]))[0] - 1,
      rtol=1e-10,
      atol=1e-10,
    )
    assert isinstance(cop.beta, float)

  # parameters_to_* at the stored parameters match the properties, and match
  # a fresh copula at other parameters.
  np.testing.assert_allclose(
    clayton.parameters_to_taildep(clayton.parameters),
    clayton.taildep,
    rtol=1e-12,
  )
  np.testing.assert_allclose(
    clayton.parameters_to_beta(clayton.parameters), clayton.beta, rtol=1e-12
  )
  other = np.array([[4.0]])
  fresh = pv.Bicop(family=pv.families.clayton, parameters=other)
  np.testing.assert_allclose(
    clayton.parameters_to_taildep(other), fresh.taildep, rtol=1e-12
  )
  np.testing.assert_allclose(
    clayton.parameters_to_beta(other), fresh.beta, rtol=1e-12
  )

  # TLL: the tail dependence coefficients are undefined (all-NaN), but
  # Blomqvist's beta is still computed from the interpolated cdf.
  rng = np.random.RandomState(5)
  u = pv.to_pseudo_obs(rng.normal(size=(200, 2)))
  tll = pv.Bicop.from_data(
    u, controls=pv.FitControlsBicop(family_set=[pv.families.tll])
  )
  assert np.isnan(tll.taildep).all()
  assert np.isfinite(tll.beta)


def test_bicop_flip_returns_swapped_copy() -> None:
  # flip() returns the argument-swapped copula and leaves the original
  # unchanged. A 90-degree Clayton is asymmetric, so the flip is a actual
  # change (rotation 90 <-> 270).
  cop = pv.Bicop.from_family(
    family=pv.families.clayton, rotation=90, parameters=np.array([[3.0]])
  )
  flipped = cop.flip()
  assert cop.rotation == 90
  assert flipped.rotation == 270
  rng = np.random.RandomState(7)
  u = rng.uniform(0.05, 0.95, size=(200, 2))
  swapped = u[:, [1, 0]]
  np.testing.assert_allclose(flipped.pdf(u), cop.pdf(swapped), rtol=1e-14)
  np.testing.assert_allclose(flipped.hfunc1(u), cop.hfunc2(swapped), rtol=1e-14)
  np.testing.assert_allclose(flipped.hfunc2(u), cop.hfunc1(swapped), rtol=1e-14)


_BICOP_SCORE_MATRIX_METHODS = ["scores", "hessian", "scores_cov"]


@pytest.mark.parametrize(
  ("family", "parameters"),
  [
    (pv.families.gaussian, np.array([[0.5]])),
    (pv.families.clayton, np.array([[2.0]])),
    (pv.families.student, np.array([[0.5], [4.0]])),
    (pv.families.bb1, np.array([[1.0], [1.5]])),
  ],
)
def test_bicop_scores_family(
  family: pv.BicopFamily, parameters: np.ndarray
) -> None:
  # Log-likelihood scores / gradient / Hessian on Bicop (vinecopulib#699),
  # mirroring the Vinecop surface for a single pair copula.
  n = 400
  cop = pv.Bicop(family=family, parameters=parameters)
  u = cop.sample(n, seeds=[1, 2, 3])
  p = cop.parameters.shape[0]

  # Shapes: scores (n, p); gradient (p,); hessian / scores_cov (p, p).
  scores = cop.scores(u)
  assert isinstance(scores, np.ndarray) and scores.shape == (n, p)
  grad = cop.gradient(u)
  assert grad.shape == (p,)
  # gradient is the observation-average of the scores.
  np.testing.assert_allclose(grad, scores.mean(axis=0), rtol=1e-10, atol=1e-12)
  assert cop.hessian(u).shape == (p, p)
  cov = cop.scores_cov(u)
  assert cov.shape == (p, p)
  np.testing.assert_allclose(cov, cov.T, rtol=1e-10, atol=1e-12)

  # hessian_full: one (p, p) matrix per observation.
  hess_full = cop.hessian_full(u)
  assert isinstance(hess_full, list) and len(hess_full) == n
  for h in hess_full:
    assert isinstance(h, np.ndarray) and h.shape == (p, p)

  # scores_full: dict carrying only the score matrix (parity with Vinecop).
  full = cop.scores_full(u)
  assert isinstance(full, dict) and set(full) == {"scores"}
  np.testing.assert_allclose(full["scores"], scores, rtol=1e-12, atol=1e-14)

  # gradient matches a central finite difference of the average log-likelihood.
  fd = np.empty(p)
  for k in range(p):
    step = 1e-6 * max(1.0, abs(cop.parameters[k, 0]))
    up, dn = cop.parameters.copy(), cop.parameters.copy()
    up[k, 0] += step
    dn[k, 0] -= step
    cop_up = pv.Bicop(family=family, parameters=up)
    cop_dn = pv.Bicop(family=family, parameters=dn)
    fd[k] = (cop_up.loglik(u) - cop_dn.loglik(u)) / (2 * step) / n
  np.testing.assert_allclose(grad, fd, rtol=1e-4, atol=1e-5)

  # Per-row parity: constant per-row parameters equal to the object's own
  # reproduce the state-based result, and threading is consistent.
  pars_const = np.tile(cop.parameters.ravel(), (n, 1))
  for method in [*_BICOP_SCORE_MATRIX_METHODS, "gradient"]:
    fn = getattr(cop, method)
    np.testing.assert_allclose(
      fn(u, parameters=pars_const), fn(u), rtol=1e-9, atol=1e-12
    )
    np.testing.assert_allclose(
      fn(u, parameters=pars_const, num_threads=2), fn(u), rtol=1e-9, atol=1e-12
    )


def test_bicop_scores_family_rejects_nonparametric_and_discrete() -> None:
  # The score family is defined for parametric, continuous pair copulas only.
  rng = np.random.RandomState(0)
  u = rng.uniform(0.05, 0.95, size=(200, 2))

  tll = pv.Bicop.from_data(
    rng.uniform(size=(300, 2)),
    controls=pv.FitControlsBicop(family_set=[pv.families.tll]),
  )
  for method in ("scores", "gradient", "hessian", "scores_cov"):
    with pytest.raises(RuntimeError):
      getattr(tll, method)(u)

  discrete = pv.Bicop(
    family=pv.families.clayton,
    parameters=np.array([[2.0]]),
    var_types=["d", "d"],
  )
  u_disc = np.column_stack([u, u])
  for method in ("scores", "gradient", "hessian", "scores_cov"):
    with pytest.raises(RuntimeError):
      getattr(discrete, method)(u_disc)


def test_bicop_family_name_and_with_var_types() -> None:
  cop = pv.Bicop.from_family(
    pv.families.gumbel, rotation=90, parameters=np.array([[2.0]])
  )
  assert cop.family_name == "Gumbel"

  discrete = pv.Bicop.from_family(
    pv.families.gaussian, parameters=np.array([[0.5]]), var_types=["d", "d"]
  )
  continuous = discrete.with_var_types()
  assert continuous.var_types == ["c", "c"]
  assert discrete.var_types == ["d", "d"]
  assert continuous.family == discrete.family
  np.testing.assert_allclose(continuous.parameters, discrete.parameters)


def test_fit_controls_bicop_selection_criterion_matches_upstream() -> None:
  # Same knob, same default, in C++, R and Python (vinecopulib#729).
  assert pv.FitControlsBicop().selection_criterion == "aic"
  assert pv.FitControlsVinecop().selection_criterion == "aic"


def test_simulate_per_row_parameters_matches_loop() -> None:
  # One parameter set per drawn observation (vinecopulib#719): row i must be
  # what a copula carrying row i's parameters would have drawn.
  n, seeds = 64, [1, 2, 3]
  cop = pv.Bicop.from_family(pv.families.clayton, parameters=np.array([[2.0]]))
  parameters = np.linspace(0.5, 6.0, n).reshape(n, 1)

  drawn = cop.sample(parameters=parameters, seeds=seeds)
  assert drawn.shape == (n, 2)

  # sample() draws an (n, 2) uniform sample and replaces the second column
  # with hinv1 of the pair; the same seeds reproduce that sample exactly.
  base = pv.utils.sample_uniform(n, 2, False, seeds)
  np.testing.assert_array_equal(drawn[:, 0], base[:, 0])

  # Row i must be what a copula carrying only row i's parameters would draw.
  for i in (0, n // 2, n - 1):
    single = pv.Bicop.from_family(
      pv.families.clayton, parameters=parameters[i : i + 1]
    )
    expected = single.hinv1(base[i : i + 1, :])
    np.testing.assert_allclose(drawn[i, 1], expected[0], rtol=1e-12)


def test_simulate_per_row_parameters_independence() -> None:
  # `indep` has no parameters, so the per-row form takes an (n, 0) matrix and
  # the row count alone fixes the sample size.
  cop = pv.Bicop.from_family(pv.families.indep)
  drawn = cop.sample(parameters=np.empty((5, 0)), seeds=[7])
  assert drawn.shape == (5, 2)


def test_simulate_requires_exactly_one_of_n_and_parameters() -> None:
  cop = pv.Bicop.from_family(pv.families.gaussian, parameters=np.array([[0.5]]))
  with pytest.raises(ValueError, match="exactly one"):
    cop.sample()
  with pytest.raises(ValueError, match="exactly one"):
    cop.sample(10, parameters=np.full((10, 1), 0.5))


def test_simulate_per_row_parameters_rejects_nonparametric() -> None:
  u = pv.to_pseudo_obs(np.random.default_rng(0).normal(size=(200, 2)))
  cop = pv.Bicop.from_data(
    u, controls=pv.FitControlsBicop(family_set=[pv.families.tll])
  )
  with pytest.raises(RuntimeError):
    cop.sample(parameters=np.full((4, 1), 0.5))


def test_simulate_per_row_parameters_must_be_two_dimensional() -> None:
  # One row per drawn observation, so a 1-d array is ambiguous: `[rho, df]` is
  # either one Student t or two malformed parameter sets.
  cop = pv.Bicop.from_family(pv.families.gaussian, parameters=np.array([[0.5]]))
  with pytest.raises(TypeError):
    cop.sample(parameters=np.array([0.1, 0.2, 0.3]))


@pytest.mark.parametrize(
  ("reorder", "fortran"),
  [(np.ascontiguousarray, False), (np.asfortranarray, True)],
)
def test_simulate_per_row_parameters_accepts_either_memory_order(
  reorder: Callable[[np.ndarray], np.ndarray], fortran: bool
) -> None:
  # `np.column_stack` and `np.stack(axis=1)` -- the obvious ways to build a
  # two-parameter array -- return C order, which must work as well as F. A
  # single-column array hides this: it is contiguous in both orders.
  n, seeds = 32, [1, 2, 3]
  parameters = reorder(
    np.column_stack([np.linspace(0.1, 0.8, n), np.linspace(3.0, 10.0, n)])
  )
  # Mutually exclusive for a multi-column array, so this pins the layout.
  assert parameters.flags.f_contiguous == fortran

  cop = pv.Bicop.from_family(
    pv.families.student, parameters=np.array([[0.5], [4.0]])
  )
  drawn = cop.sample(parameters=parameters, seeds=seeds)
  assert drawn.shape == (n, 2)

  reference = cop.sample(parameters=np.asfortranarray(parameters), seeds=seeds)
  np.testing.assert_array_equal(drawn, reference)


def test_simulate_positional_signature_is_unchanged() -> None:
  # `sample(n, qrng, seeds)` predates the per-row overload and must keep
  # meaning what it meant.
  cop = pv.Bicop.from_family(pv.families.gaussian, parameters=np.array([[0.5]]))
  assert cop.sample(12, False, seeds=[1, 2]).shape == (12, 2)


def _discrete_pair(
  n: int = 500, seed: int = 3
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  """A (continuous, discrete) pair as F(x) and the discrete left limit."""
  rng = np.random.default_rng(seed)
  z = rng.normal(size=(n, 2))
  u1 = pv.to_pseudo_obs(z[:, [0]])[:, 0]
  counts = np.floor(np.abs(z[:, 1]) * 3).astype(int)
  hi = (counts + 1) / (counts.max() + 2)
  lo = counts / (counts.max() + 2)
  return u1, hi, lo


def test_from_data_accepts_the_compact_discrete_layout() -> None:
  # `from_data` took a statically two-column matrix, so a discrete pair — which
  # needs a left-limit column — could not be passed at all (vinecopulib#729).
  u1, hi, lo = _discrete_pair()
  compact = np.column_stack([u1, hi, lo])
  cop = pv.Bicop.from_data(compact, var_types=["c", "d"])
  assert cop.var_types == ["c", "d"]


def test_from_data_discrete_layouts_agree() -> None:
  # Expanded (n, 4) and compact (n, 2 + k) describe the same data.
  u1, hi, lo = _discrete_pair()
  compact = np.column_stack([u1, hi, lo])
  expanded = np.column_stack([u1, hi, u1, lo])

  controls = pv.FitControlsBicop(family_set=[pv.families.gaussian])
  a = pv.Bicop.from_data(compact, controls=controls, var_types=["c", "d"])
  b = pv.Bicop.from_data(expanded, controls=controls, var_types=["c", "d"])
  np.testing.assert_allclose(a.parameters, b.parameters, rtol=1e-12)


def test_from_data_still_accepts_plain_continuous_input() -> None:
  # Widening the accepted type must not change how ordinary input converts.
  # nanobind's ndarray caster takes arrays, not sequences, on every
  # Eigen-typed argument in this package; `from_data` is no exception.
  rng = np.random.default_rng(0)
  u = pv.to_pseudo_obs(rng.normal(size=(200, 2)))
  assert isinstance(pv.Bicop.from_data(u), pv.Bicop)
  assert isinstance(pv.Bicop.from_data(np.asarray(u, order="C")), pv.Bicop)
  with pytest.raises(TypeError):
    pv.Bicop.from_data(u.tolist())


def test_loglik_rejects_a_wrong_column_count() -> None:
  # vinecopulib#729 made `loglik` validate the column count instead of
  # reading past the data.
  rng = np.random.default_rng(1)
  u = pv.to_pseudo_obs(rng.normal(size=(100, 2)))
  cop = pv.Bicop.from_data(u)
  with pytest.raises(RuntimeError):
    cop.loglik(np.column_stack([u, u[:, [0]], u[:, [1]], u[:, [0]]]))


@pytest.mark.parametrize(
  "family_set",
  [
    [pv.families.gaussian],
    [pv.families.clayton, pv.families.gumbel, pv.families.frank],
    [pv.families.tll],
  ],
)
def test_a_weighted_fit_is_flip_equivariant(
  family_set: list[pv.BicopFamily],
) -> None:
  """Weighting the rows cannot depend on which argument comes first.

  A weighted fit on ``(u1, u2)`` must equal the weighted fit on ``(u2, u1)`` up
  to a flip: the weights are per-observation, so exchanging the two arguments
  exchanges the fitted copula's arguments and nothing else. It is the pair-copula
  analog of the argument-order invariance ``find_latent_sample`` guarantees,
  and nothing pinned it.
  """
  rng = np.random.default_rng(0)
  u = pv.to_pseudo_obs(
    rng.multivariate_normal([0.0, 0.0], [[1.0, 0.6], [0.6, 1.0]], 600)
  )
  controls = pv.FitControlsBicop(family_set=family_set, num_threads=1)
  controls.weights = rng.uniform(0.2, 2.0, u.shape[0])

  forward = pv.Bicop.from_data(u, controls=controls)
  reversed_ = pv.Bicop.from_data(u[:, [1, 0]], controls=controls)
  flipped = forward.flip()

  assert flipped.family == reversed_.family
  assert flipped.rotation == reversed_.rotation
  np.testing.assert_allclose(
    np.asarray(flipped.parameters),
    np.asarray(reversed_.parameters),
    rtol=1e-10,
    atol=1e-10,
  )
  # And the densities agree with the arguments exchanged.
  grid = rng.uniform(0.02, 0.98, size=(200, 2))
  np.testing.assert_allclose(
    forward.pdf(grid), reversed_.pdf(grid[:, [1, 0]]), rtol=1e-10, atol=1e-10
  )
