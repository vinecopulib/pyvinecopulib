import numpy as np
import pytest

import pyvinecopulib as pv


def test_sample_uniform() -> None:
  assert isinstance(pv.utils.sample_uniform(10, 2, False, [1, 2]), np.ndarray)


@pytest.mark.parametrize("generator", [pv.utils.ghalton, pv.utils.sobol])
@pytest.mark.parametrize(("n", "d"), [(0, 2), (2, 0), (0, 0)])
def test_qrng_rejects_empty_shapes(generator, n: int, d: int) -> None:
  # Before vinecopulib#719 only the pseudo-random branch validated, so the
  # quasi-random generators wrote out of bounds for a zero extent.
  with pytest.raises(RuntimeError):
    generator(n, d)
