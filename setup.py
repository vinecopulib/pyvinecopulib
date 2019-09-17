import os
import re
import sys
import platform
import subprocess
from glob import glob

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


if not os.path.isdir("lib/boost"):
    import tarfile
    boost_archive = glob(os.path.join("lib", "boost*.tar.gz"))[0]
    tar = tarfile.open(boost_archive)
    tar.extractall(path='lib')
    tar.close()


def get_sources(paths):
    sources = []
    for path in paths:
        [sources.append(y) for x in os.walk(path)
         for y in glob(os.path.join(x[0], '*'))
         if not os.path.isdir(y)]
    return sources


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        paths = ['boost', 'eigen/unsupported', 'eigen/Eigen', 'pybind11',
                 'vinecopulib/include', 'wdm/include']
        sources = get_sources(['lib/' + path for path in paths])
        Extension.__init__(self, name, sources=sources)
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='pyvinecopulib',
    version='0.0.1',
    install_requires=requirements,
    author='Thomas Nagler and Thibault Vatter',
    author_email='info@vinecopulib.org',
    description='A python interface to vinecopulib',
    long_description='TODO',
    url="https://github.com/pyvinecopulib/",
    ext_modules=[CMakeExtension('pyvinecopulib')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
