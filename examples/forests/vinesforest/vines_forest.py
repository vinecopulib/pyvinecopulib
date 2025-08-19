from typing import Optional

import numpy as np
from sklearn.model_selection import train_test_split

import pyvinecopulib as pv

ALGORITHMS = ["jck", "wilson"]


class VinesForest:
  def __init__(
    self,
    n_vines: int = 100,
    train_val_split: float = 0.8,
    algorithm: str = "jck",
    controls: pv.FitControlsVinecop = pv.FitControlsVinecop(
      family_set=[pv.tll], num_threads=1
    ),
    seed: int = 42,
  ) -> None:
    self.n_vines = n_vines
    self.algorithm = algorithm
    self.train_val_split = train_val_split
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
    self.seed = seed

  def fit(
    self,
    data: np.ndarray,
    data_val: Optional[np.ndarray] = None,
  ) -> None:
    d = data.shape[1]

    # Set seed and create random number generator
    rng = np.random.default_rng(self.seed)

    # For the results
    self.vines = []
    logliks_forest = []

    # Save the initial control seed
    initial_controls_seed = self.controls.seeds

    # If no validation data is provided, we split the data into train and validation
    if data_val is None:
      train, val = train_test_split(
        data, train_size=self.train_val_split, random_state=self.seed
      )
    else:
      if data_val.shape[1] != d:
        raise ValueError(
          f"Validation data must have the same number of dimensions as training data. Expected {d}, got {data_val.shape[1]}"
        )
      # If validation data is provided, we use it directly
      train = data
      val = data_val
    self.vine_dissman = pv.Vinecop.from_data(data=train, controls=self.controls)
    self.loglik_dissman = np.log(self.vine_dissman.pdf(val))

    for _ in range(self.n_vines):
      # Seeding
      seeds = [int(rng.integers(0, 2**31 - 1)) for _ in range(2)]
      self.controls.seeds = seeds

      # If the algorithm is jck, we simulate the structure
      structure = (
        pv.RVineStructure.simulate(d, seeds=seeds)
        if self.algorithm == "jck"
        else None
      )

      # Fit the vine
      vine = pv.Vinecop.from_data(
        data=train, structure=structure, controls=self.controls
      )

      # Scoring per observation
      loglik_vine = np.log(vine.pdf(val))

      # Saving the vine and its per observation validation loglikelihood
      self.vines.append(vine)
      logliks_forest.append(loglik_vine)

    # Store the loglikelihoods of the forest
    self.logliks_forest = np.array(logliks_forest)

    # Restore the initial control seed
    self.controls.seeds = initial_controls_seed

  # def loglik(
  #   self, data: Optional[np.ndarray] = None, per_observation: bool = False
  # ) -> np.ndarray:
  #   if not self.vines:
  #     raise ValueError("VinesForest not fitted yet")
  #   if data is None:
  #     data = np.ndarray(shape=(0, 0), dtype=np.float64)
  #   if not per_observation:
  #     return np.array([vine.loglik(data) for vine in self.vines])
  #   else:
  #     return np.log(np.concatenate([vine.pdf(data) for vine in self.vines]))
