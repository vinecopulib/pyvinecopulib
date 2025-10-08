import uuid
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

# # Per-test unique directory, xdist-friendly (each worker gets its own base dir)
# @pytest.fixture(scope="function")
# def test_dump_folder(
#   tmp_path_factory: pytest.TempPathFactory, request: pytest.FixtureRequest
# ) -> str:
#   worker = os.getenv("PYTEST_XDIST_WORKER", "worker0")
#   base = tmp_path_factory.mktemp(f"dump-{worker}")
#   # also isolate by test node to avoid clashes inside same worker
#   d = base / request.node.name.replace(os.sep, "_")
#   d.mkdir(parents=True, exist_ok=True)
#   return str(os.fspath(d))


# If you prefer a unique file path directly:
@pytest.fixture
def unique_json_path(tmp_path: Path, request: pytest.FixtureRequest) -> Path:
  return tmp_path / f"{request.node.name}-{uuid.uuid4().hex}.json"
