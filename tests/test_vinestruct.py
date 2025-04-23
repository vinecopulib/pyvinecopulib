import pyvinecopulib as pv

from .helpers import compare_rvinestructure


def test_rvinestructure(test_dump_folder: str) -> None:
  # Test RVineStructure class
  rvine = pv.RVineStructure(5)
  assert isinstance(rvine, pv.RVineStructure)
  assert rvine.dim == 5
  assert rvine.matrix.shape == (5, 5)
  assert rvine.trunc_lvl == 4
  assert rvine.order == list(range(1, 6))

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
