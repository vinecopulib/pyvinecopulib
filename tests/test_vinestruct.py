import pyvinecopulib as pv

from .helpers import compare_rvinestructure


def test_rvinestructure(test_dump_folder: str) -> None:
  d = 5
  # Test RVineStructure class
  rvine = pv.RVineStructure(d)
  assert isinstance(rvine, pv.RVineStructure)
  assert rvine.dim == d
  assert rvine.trunc_lvl == d - 1
  assert rvine.order == list(range(1, d + 1))
  matrix = rvine.matrix
  assert isinstance(matrix, np.ndarray)
  assert matrix.shape == (d, d)
  assert matrix.dtype == np.uint64
  assert np.all(np.logical_and(matrix >= 0, matrix <= d))

  # Should be the same as the previous one
  dvine = pv.DVineStructure(rvine.order)
  compare_rvinestructure(dvine, rvine, True)

  # Test to_json and from_json
  new_rvine = pv.RVineStructure.from_json(rvine.to_json())
  compare_rvinestructure(rvine, new_rvine)
  filename = test_dump_folder + "/test_rvine.json"
  rvine.to_file(filename)
  new_rvine = pv.RVineStructure.from_file(filename)
  compare_rvinestructure(rvine, new_rvine)
