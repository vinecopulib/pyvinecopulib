import json
import os

import numpy as np
import pytest

import pyvinecopulib as pv


def test_bicop(unique_json_path: str) -> None:
  bicop = pv.Bicop()

  # Test default initialization
  assert bicop.family == pv.families.indep
  assert bicop.rotation == 0
  assert bicop.parameters.shape == (0, 0)
  assert bicop.var_types == ["c", "c"]

  # Test initialization with arguments
  data = np.array([[0.1, 0.2], [0.3, 0.4]])
  controls = pv.FitControlsBicop()
  bicop = pv.Bicop.from_data(data, controls)

  assert bicop.family == pv.families.indep
  assert bicop.rotation == 0
  assert bicop.parameters.shape == (0, 0)
  assert bicop.var_types == ["c", "c"]

  # Test to_json method
  new_bicop = pv.Bicop.from_json(bicop.to_json())
  assert bicop.family == new_bicop.family
  assert bicop.rotation == new_bicop.rotation
  assert bicop.parameters.shape == new_bicop.parameters.shape
  assert bicop.var_types == new_bicop.var_types
  filename = os.fspath(unique_json_path)
  bicop.to_file(filename)
  assert os.path.exists(filename)
  new_bicop = pv.Bicop.from_file(filename)
  assert bicop.family == new_bicop.family
  assert bicop.rotation == new_bicop.rotation
  assert bicop.parameters.shape == new_bicop.parameters.shape
  assert bicop.var_types == new_bicop.var_types

  # A non-.cbor filename keeps writing JSON text (backwards compatibility) ...
  with open(filename, encoding="utf-8") as f:
    json.load(f)

  # ... while a ``.cbor`` filename selects binary CBOR (vinecopulib#684).
  cbor_filename = filename.removesuffix(".json") + ".cbor"
  bicop.to_file(cbor_filename)
  new_bicop = pv.Bicop.from_file(cbor_filename)
  assert bicop.family == new_bicop.family
  assert bicop.rotation == new_bicop.rotation
  assert bicop.parameters.shape == new_bicop.parameters.shape
  assert bicop.var_types == new_bicop.var_types
  with open(cbor_filename, "rb") as f:
    assert f.read(1) != b"{"

  # Test properties
  bicop = pv.Bicop(family=pv.families.gumbel, rotation=90)
  bicop.rotation = 0
  assert bicop.rotation == 0
  with pytest.raises(RuntimeError):
    bicop.rotation = 45

  bicop.parameters = np.array([[3.0]])
  assert bicop.parameters.shape == (1, 1)
  assert bicop.parameters[0, 0] == 3.0

  bicop.var_types = ["d", "d"]
  assert bicop.var_types == ["d", "d"]

  # Test read-only properties
  assert isinstance(bicop.tau, float)
  assert bicop.npars == 1
  with pytest.raises(AttributeError):
    setattr(bicop, "npars", 2)

  # Test passing a single row of data (#169 & #170 fix)
  bicop.var_types = ["c", "c"]
  u = np.array([[0.1, 0.2]])
  d = bicop.pdf(u)
  assert isinstance(d, np.ndarray) and d.shape == (1,)

  # Test loglik method
  u = np.array([[0.1, 0.2], [0.3, 0.4]])
  loglik = bicop.loglik(u)
  assert isinstance(loglik, float)

  # Test aic method
  aic = bicop.aic(u)
  assert isinstance(aic, float)

  # Test bic method
  bic = bicop.bic(u)
  assert isinstance(bic, float)

  # Test mbic method
  psi0 = 0.9
  mbic = bicop.mbic(u, psi0)
  assert isinstance(mbic, float)

  # Test __repr__ method
  assert isinstance(repr(bicop), str)

  # Test str method
  assert isinstance(str(bicop), str)

  # Test parameters_to_tau method
  parameters = np.array([[0.5, 0.6], [0.7, 0.8]])
  tau = bicop.parameters_to_tau(parameters)
  assert isinstance(tau, float)

  # Test tau_to_parameters method
  tau = 0.5
  parameters = bicop.tau_to_parameters(tau)
  assert isinstance(parameters, np.ndarray)

  # Test parameters_lower_bounds method
  lower_bounds = bicop.parameters_lower_bounds
  assert isinstance(lower_bounds, np.ndarray)
  assert lower_bounds == np.array([1.0])

  # Test parameters_upper_bounds method
  upper_bounds = bicop.parameters_upper_bounds
  assert isinstance(upper_bounds, np.ndarray)
  assert upper_bounds == np.array([50.0])

  for method in ["pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"]:
    values = getattr(bicop, method)(u)
    assert isinstance(values, np.ndarray)
    assert values.shape == (2,)

  # Test simulate method
  n = 100
  qrng = False
  seeds: list[int] = []
  samples = bicop.simulate(n, qrng, seeds)
  assert samples.shape == (n, 2)

  # Test fit method
  controls = pv.FitControlsBicop()
  bicop.fit(u, controls)

  # Test select method
  controls = pv.FitControlsBicop()
  bicop.select(u, controls)


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
    values = fn(u, pars)
    assert isinstance(values, np.ndarray)
    assert values.shape == (n,)

    # Identity parity: constant per-row parameters equal to the object's own
    # parameters reproduce the state-based (single-argument) result.
    np.testing.assert_allclose(fn(u, pars_const), fn(u), rtol=1e-9, atol=1e-12)

    # Threading must not change the result.
    np.testing.assert_allclose(
      fn(u, pars, num_threads=2), values, rtol=1e-12, atol=1e-14
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
  ll = cop.loglik(u, pars)
  assert isinstance(ll, float)
  np.testing.assert_allclose(ll, np.nansum(np.log(cop.pdf(u, pars))), rtol=1e-9)
  np.testing.assert_allclose(
    cop.loglik(u, pars_const), cop.loglik(u), rtol=1e-9, atol=1e-12
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
    tll.pdf(u, np.ones((n, 1)))

  # Wrong number of parameter rows (must equal u.rows()).
  with pytest.raises(RuntimeError):
    cop.pdf(u, pars[:-1])

  # Wrong number of parameter columns (must equal the family's parameter count).
  with pytest.raises(RuntimeError):
    cop.pdf(u, np.ones((n, 2)))

  # Non-finite parameters are rejected.
  pars_nan = pars.copy()
  pars_nan[0, 0] = np.nan
  with pytest.raises(RuntimeError):
    cop.pdf(u, pars_nan)

  # Out-of-bounds parameters are rejected.
  pars_oob = pars.copy()
  pars_oob[0, 0] = -5.0
  with pytest.raises(RuntimeError):
    cop.pdf(u, pars_oob)


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
    values = fn(u, sel, pars)
    assert isinstance(values, np.ndarray) and values.shape == (n,)

    # Constant per-row parameters reproduce the state-based call.
    np.testing.assert_allclose(
      fn(u, sel, pars_const), fn(u, sel), rtol=1e-9, atol=1e-12
    )

    # Threading must not change the result.
    np.testing.assert_allclose(
      fn(u, sel, pars, num_threads=2), values, rtol=1e-12, atol=1e-14
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

  # Parameter selectors, second-order derivatives, and genuine h-function
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
  # deliberately.
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
