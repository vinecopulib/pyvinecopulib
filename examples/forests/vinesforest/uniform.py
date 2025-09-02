import numpy as np


class Uniform:
  def __init__(self, xmin=0, xmax=1, type="continuous"):
    if type != "continuous":
      raise ValueError("Uniform distribution only supports continuous type.")
    self.xmin = xmin
    self.xmax = xmax
    self.type = type

  def fit(self, x):
    pass

  def cdf(self, x):
    return np.clip((x - self.xmin) / (self.xmax - self.xmin), 0, 1)

  def pdf(self, u):
    return np.where(
      (u >= self.xmin) & (u <= self.xmax), 1 / (self.xmax - self.xmin), 0
    )
