# Import the required libraries
import numpy as np
import pyvinecopulib as pv
from scipy.stats import poisson

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
fit_cont = pv.Vinecop(u, controls=controls)
str(fit_cont)

# Model for discrete data
# Transform to Poisson margins
x_poisson = poisson.ppf(u, 1)

# using Poisson(1) transformation
u_disc = np.hstack((poisson.cdf(x_poisson, 1), poisson.cdf(x_poisson - 1, 1)))

# Fit vine copula model for discrete data
fit_disc = pv.Vinecop(u_disc, var_types=["d"] * 5, controls=controls)
str(fit_disc)

# Model for mixed data
# Transform first variable to Poisson margin
x_poisson_mixed = poisson.ppf(u[:, 0], 1)
u_disc_mixed = np.hstack(
  (
    poisson.cdf(x_poisson_mixed, 1).reshape(-1, 1),
    u[:, 1:5],
    poisson.cdf(x_poisson_mixed - 1, 1).reshape(-1, 1),
  )
)

# Fit vine copula model for mixed data
fit_mixed = pv.Vinecop(
  u_disc_mixed, var_types=["d"] + ["c"] * 4, controls=controls
)
str(fit_mixed)
