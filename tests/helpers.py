import numpy as np


def random_data(d=5, n=1000):
  # Simulate some data
  np.random.seed(1234)  # seed for the random generator
  mean = np.random.normal(size=d)  # mean vector
  cov = np.random.normal(size=(d, d))  # covariance matrix
  cov = np.dot(cov.transpose(), cov)  # make it non-negative definite
  x = np.random.multivariate_normal(mean, cov, n)
  return x


def compare_properties(obj1, obj2, attrs, subclass=False):
  if subclass:
    assert issubclass(
      type(obj1), type(obj2)
    ), "Objects must be of the same type"
  else:
    assert type(obj1) is type(obj2), "Objects must be of the same type"
  for attr in attrs:
    val1 = getattr(obj1, attr)
    val2 = getattr(obj2, attr)
    if isinstance(val1, np.ndarray):
      assert isinstance(val2, np.ndarray) and np.array_equal(
        val1, val2
      ), f"Mismatch in {attr}: {val1} != {val2}"
    else:
      assert val1 == val2, f"Mismatch in {attr}: {val1} != {val2}"


def compare_bicop(cop1, cop2):
  attrs = ["family", "rotation", "parameters", "var_types"]
  compare_properties(cop1, cop2, attrs)


def compare_vinecop(cop1, cop2):
  attrs = ["dim", "trunc_lvl", "order", "matrix"]
  compare_properties(cop1, cop2, attrs)

  d = cop1.dim
  trunc_lvl = cop1.trunc_lvl
  for t in range(trunc_lvl):
    for e in range(d - t - 1):
      bicop1 = cop1.get_pair_copula(t, e)
      bicop2 = cop2.get_pair_copula(t, e)
      compare_bicop(bicop1, bicop2)


def compare_rvinestructure(struct1, struct2, subclass=False):
  attrs = ["dim", "order", "trunc_lvl", "matrix"]
  compare_properties(struct1, struct2, attrs, subclass)
