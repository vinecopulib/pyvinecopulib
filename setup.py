import os
import subprocess
from pathlib import Path
from sys import platform

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup


def find_include_dirs():
  if platform == "linux" or platform == "linux2" or platform == "darwin":
    try:
      # Run CMake and capture its output
      cmake_output = subprocess.check_output(
        ["cmake", "."], stderr=subprocess.STDOUT, text=True
      )

      # Extract the include directories from CMake output
      for line in cmake_output.splitlines():
        if line.startswith("-- EIGEN3_INCLUDE_DIRS="):
          eigen_include = line.split("=")[1].strip()
        elif line.startswith("-- Boost_INCLUDE_DIRS="):
          boost_include = line.split("=")[1].strip()

    except subprocess.CalledProcessError as e:
      raise RuntimeError(
        f"Could not find Boost or Eigen using CMake. Ensure they are installed. Error: {e.output}"
      )

  elif platform == "win32":
    boost_include = os.environ.get("BOOST_INCLUDE_DIR")
    eigen_include = os.environ.get("EIGEN_INCLUDE_DIR")

    if not boost_include or not eigen_include:
      raise RuntimeError(
        "BOOST_INCLUDE_DIR and EIGEN_INCLUDE_DIR environment variables must be set."
      )

  include_dirs = [
    boost_include,
    eigen_include,
  ]  # , os.path.join(eigen_include, "unsupported")]

  return include_dirs


include_dirs = find_include_dirs()
include_dirs.append("lib/vinecopulib/include")
include_dirs.append("lib/wdm/include")

setup(
  name="pyvinecopulib",
  long_description=(Path("setup.py").parent / "README.md").read_text(),
  long_description_content_type="text/markdown",
  ext_modules=[
    Pybind11Extension(
      "pyvinecopulib",
      ["src/main.cpp"],
      include_dirs=include_dirs,
      language="c++",
      cxx_std=17,
    )
  ],
  cmdclass={"build_ext": build_ext},
  zip_safe=False,
)
