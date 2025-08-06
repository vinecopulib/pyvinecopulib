from typing import Optional

import numpy as np

import pyvinecopulib as pv

ALGORITHMS = ["jck", "wilson"]


class VinesForest:
  def __init__(
    self,
    n_vines: int = 100,
    algorithm: str = "jck",
    controls: pv.FitControlsVinecop = pv.FitControlsVinecop(
      family_set=[pv.tll], num_threads=1
    ),
  ) -> None:
    self.n_vines = n_vines
    self.algorithm = algorithm
    if algorithm not in ALGORITHMS:
      raise ValueError(f"Algorithm must be one of {ALGORITHMS}")
    if algorithm == "wilson":
      controls.tree_algorithm = "random_weighted"
    self.controls = controls
    # if self.controls.num_threads != 1:
    #   self.controls.num_threads = 1
    #   warnings.warn(
    #     "num_threads > 1 is not supported by VinesForest. Setting num_threads=1."
    #   )
    self.vines: list[tuple[int, pv.Vinecop]] = []

  def fit(
    self,
    data: np.ndarray,
    trunc_lvl: int = int(2**31 - 1),
    seed: int = 42,
  ) -> None:
    # Adjust and set trunc_lvl
    dim = data.shape[1]
    self.controls.trunc_lvl = min(dim, trunc_lvl)

    # Set seed and create random number generator
    rng = np.random.default_rng(seed)

    # For the results
    self.vines = []

    # These seeds ensure that the sampled vines are different
    seeds = [int(rng.integers(0, 2**31 - 1)) for _ in range(self.n_vines)]

    for seed_i in seeds:
      rng = np.random.default_rng(seed_i)
      seeds_i = [int(rng.integers(0, 2**31 - 1)) for _ in range(10)]
      if self.algorithm == "jck":
        structure = pv.RVineStructure.simulate(data.shape[1], seeds=seeds_i)
        vine = pv.Vinecop.from_data(
          data=data, structure=structure, controls=self.controls
        )
      else:
        self.controls.seeds = seeds_i
        vine = pv.Vinecop.from_data(data=data, controls=self.controls)
      self.vines.append((seed_i, vine))

  def loglik(
    self, data: Optional[np.ndarray] = None, per_observation: bool = False
  ) -> np.ndarray:
    if not self.vines:
      raise ValueError("VinesForest not fitted yet")
    if data is None:
      data = np.ndarray(shape=(0, 0), dtype=np.float64)
    if not per_observation:
      return np.array([vine.loglik(data) for _, vine in self.vines])
    else:
      return np.concatenate([vine.pdf(data) for _, vine in self.vines]).log()
