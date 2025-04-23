import numpy as np
from numpy.typing import NDArray

import pyvinecopulib as pv


def test_simulate_uniform() -> None:
  assert isinstance(
    pv.simulate_uniform(10, 2, False, [1, 2]), NDArray[np.float64]
  )
