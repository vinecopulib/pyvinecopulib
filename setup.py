import os
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

print("Present working directory: ", os.getcwd())
print("Environment variables: ", os.environ)


def find_include_dirs():
  print("Environment variables: ", os.environ)
  boost_include = os.environ.get("Boost_INCLUDE_DIR")
  if not boost_include:
    raise RuntimeError("Boost_INCLUDE_DIR environment variables must be set.")

  eigen_include = os.environ.get("EIGEN3_INCLUDE_DIR")
  if not eigen_include:
    raise RuntimeError("EIGEN3_INCLUDE_DIR environment variables must be set.")

  # Check if directories exist and are accessible
  if not os.path.isdir(boost_include):
    raise RuntimeError(f"Boost include directory not found: {boost_include}")
  if not os.path.isdir(eigen_include):
    raise RuntimeError(f"Eigen include directory not found: {eigen_include}")

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
