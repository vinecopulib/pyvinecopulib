import pyvinecopulib as pv
import numpy as np

#  def test_version():
#  assert pv.__version__ == '0.0.1'


def test_simulate_uniform():
    assert isinstance(
        pv.simulate_uniform(10, 2, False, [1, 2]),
        np.ndarray
    )
