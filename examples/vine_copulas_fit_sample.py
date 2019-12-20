# Import the required libraries
import pyvinecopulib as pv
import numpy as np

# Simulate some data
np.random.seed(1234)  # seed for the random generator
n = 1000  # number of observations
d = 5  # the dimension
mean = np.random.normal(size=d)  # mean vector
cov = np.random.normal(size=(d, d))  # covariance matrix
cov = np.dot(cov.transpose(), cov)  # make it non-negative definite
x = np.random.multivariate_normal(mean, cov, n)

# Transform copula data using the empirical distribution
u = pv.to_pseudo_obs(x)

# Fit a Gaussian vine
# (i.e., properly specified since the data is multivariate normal)
controls = pv.FitControlsVinecop(family_set=[pv.BicopFamily.gaussian])
cop = pv.Vinecop(u, controls=controls)

# Sample from the copula
n_sim = 1000
u_sim = cop.simulate(n_sim, seeds=[1, 2, 3, 4])

# Transform back simulations to the original scale
x_sim = np.asarray([np.quantile(x[:, i], u_sim[:, i]) for i in range(0, d)])

# Both the mean and covariance matrix look ok!
[mean, np.mean(x_sim, 1)]
[cov, np.cov(x_sim)]
