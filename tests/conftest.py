import os
import shutil
from typing import Iterator

import pytest


@pytest.fixture
def test_dump_folder() -> Iterator[str]:
  folder = "test_dump"
  os.makedirs(folder, exist_ok=True)
  yield folder
  shutil.rmtree(folder, ignore_errors=True)
