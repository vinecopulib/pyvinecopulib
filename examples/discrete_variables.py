# Import the required libraries
import math

import numpy as np

import pyvinecopulib as pv

# Simulate dummy data
np.random.seed(1234)  # seed for the random generator
x = np.random.normal(size=(30, 1)) * np.ones((30, 5)) + 0.5 * np.random.normal(
  size=(30, 5)
)

# Convert to pseudo-observations
u = pv.to_pseudo_obs(x)

# Some fit controls
controls = pv.FitControlsVinecop(family_set=[pv.BicopFamily.gaussian])


# A continuous example
fit_cont = pv.Vinecop.from_data(u, controls=controls)
str(fit_cont)


# Model for discrete data
# Transform to Poisson margins
# Percent Point Function (Inverse CDF, PPF)
@np.vectorize
def poisson_ppf(p, mu):
  cumulative_prob = 0.0
  k = 0
  while cumulative_prob < p:
    cumulative_prob += math.exp(-mu) * (mu**k) / math.factorial(k)
    if cumulative_prob >= p:
      return k
    k += 1
  return k  # In case p is exactly 1


x_poisson = poisson_ppf(u, 1)


# Using Poisson(1) transformation
# Cumulative Distribution Function (CDF)
@np.vectorize
def poisson_cdf(k, mu):
  return sum(
    math.exp(-mu) * (mu**i) / math.factorial(i) for i in range(int(k) + 1)
  )


u_disc = np.hstack((poisson_cdf(x_poisson, 1), poisson_cdf(x_poisson - 1, 1)))

# Fit vine copula model for discrete data
fit_disc = pv.Vinecop.from_data(u_disc, var_types=["d"] * 5, controls=controls)
str(fit_disc)

# Model for mixed data
# Transform first variable to Poisson margin
x_poisson_mixed = poisson_ppf(u[:, 0], 1)
u_disc_mixed = np.hstack(
  (
    poisson_cdf(x_poisson_mixed, 1).reshape(-1, 1),
    u[:, 1:5],
    poisson_cdf(x_poisson_mixed - 1, 1).reshape(-1, 1),
  )
)

# Fit vine copula model for mixed data
fit_mixed = pv.Vinecop.from_data(
  u_disc_mixed, var_types=["d"] + ["c"] * 4, controls=controls
)
str(fit_mixed)
