import numpy as np
import pyvinecopulib as pv


def test_simulate_uniform():
    assert isinstance(
        pv.simulate_uniform(10, 2, False, [1, 2]),
        np.ndarray
    )
