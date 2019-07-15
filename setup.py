# /*============================================================================
#
#  PYVINECOPULIB: A python interface to vinecopulib.
#
#  Copyright (c) University College London (UCL). All rights reserved.
#
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
#  See LICENSE.txt in the top level directory for details.
#
# ============================================================================*/

# Note: The whole premise of this setup.py is that this project is primarily
# a C++ project, for C++ developers. So, all CMake-ing, and Building is
# done in a pre-existing and pre-configured C++ build directory.
# The only CMake variables that this script sets are:
#   PYVINECOPULIB_PYTHON_MODULE_NAME to give the output python module the correct name.
#   PYVINECOPULIB_PYTHON_OUTPUT_DIRECTORY to produce the output python module in the correct place
# Furthermore, the Python build environment must not pick up a different version of CMake.
# The required CMake version is set in CMakeLists.txt as a C++ developer would expect.

import os
import platform
import re
import subprocess
import six
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
import versioneer

# Get the long description, as top-level README.md
with open('README.md') as f:
    long_description = f.read()

# Get the top-level folder name of this project.
dir_path = os.path.dirname(os.path.abspath(__file__))
dir_name = os.path.split(dir_path)[1]


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.source_dir = os.path.abspath(sourcedir)
        self.super_build_dir = os.path.join(self.source_dir, 'build')
        self.build_dir = os.path.join(self.super_build_dir, 'PYVINECOPULIB-build')


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
        six.print_("CMake version in python build:" + str(cmake_version))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        ext_dir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(self.distribution.get_name())))
        cmake_args = ['-DPYVINECOPULIB_PYTHON_OUTPUT_DIRECTORY:PATH=' + ext_dir,
                      '-DPYVINECOPULIB_PYTHON_MODULE_NAME:STRING=' + ext.name
                      ]
        build_args = []

        six.print_("build_extension:name=" + str(ext.name))
        six.print_("build_extension:ext_dir=" + str(ext_dir))
        six.print_("self.distribution.get_name()=" + str(self.distribution.get_name()))

        cfg = 'Debug' if self.debug else 'Release'

        if platform.system() == "Windows":
            cmake_args += ['--config ' + cfg]
            build_args += ['--config', cfg]
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        if os.environ.get('COMPILER') is not None:
            cmake_args += ['-G', str(os.environ.get('COMPILER'))]

        env = os.environ.copy()

        subprocess.check_call(['cmake'] + cmake_args + [ext.source_dir], cwd=ext.build_dir, env=env)
        subprocess.check_call(['cmake'] + ['--build', '.'] + build_args, cwd=ext.build_dir)


setup(
    # Must match python module name in your c++ code, or else you end
    # up with two dynamically linked libraries inside one wheel.
    name=('pyvinecopulib'),

    # Must match the version number in CMakeLists.txt.
    # We could try to parse the CMakeLists.txt file, but lets keep it simple.
    version=versioneer.get_version(),
    author='Matt Clarkson',
    author_email='m.clarkson@ucl.ac.uk',
    description='A template project, to enable people to build nicely structured C++ projects.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    ext_modules=[CMakeExtension('', sourcedir=dir_path)],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    license='BSD-3 license',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Software Development',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
    ],

    keywords='C++ cmake catch project template',

    install_requires=[
        'six>=1.10',
        'numpy>=1.11',
    ],
)
