import os
import tarfile
from glob import glob
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# Extract boost if not already done.
if not os.path.isdir('lib/boost'):
  boost_archive = glob(os.path.join('lib', 'boost*.tar.gz'))[0]
  tar = tarfile.open(boost_archive)
  tar.extractall(path='lib/boost')
  tar.close()

setup(
    name='pyvinecopulib',
    long_description=(Path('setup.py').parent / 'README.md').read_text(),
    long_description_content_type='text/markdown',
    ext_modules=[
        Pybind11Extension('pyvinecopulib', ['src/main.cpp'],
                          include_dirs=[
                              'lib/boost', 'lib/eigen',
                              'lib/eigen/unsupported',
                              'lib/vinecopulib/include', 'lib/wdm/include'
                          ],
                          language='c++',
                          cxx_std=17)
    ],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
