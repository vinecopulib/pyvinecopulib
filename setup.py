import os
import sys
from glob import glob

import setuptools
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as file:
        file.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([file.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False

    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag."""
    flags = ['-std=c++14', '-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag):
            return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])

        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' %
                        self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))

            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' %
                        self.distribution.get_version())

        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        build_ext.build_extensions(self)


def extract_boost():
    """Extract boost if not already done."""

    if not os.path.isdir("lib/boost"):
        import tarfile
        boost_archive = glob(os.path.join("lib", "boost*.tar.gz"))[0]
        tar = tarfile.open(boost_archive)
        tar.extractall(path='lib/boost')
        tar.close()


def get_requirements():
    """Read and return requirements."""
    with open('requirements.txt') as file:
        requirements = file.read().splitlines()

    return requirements


def get_long_description():
    """Read and return the long description."""
    with open("README.md", "r") as file:
        long_description = file.read()

    return long_description


def get_include_paths():
    """Return the long description."""
    include_dirs = ['boost', 'eigen', 'eigen/unsupported',
                    'pybind11/include', 'vinecopulib/include', 'wdm/include']

    return ['lib/' + path for path in include_dirs]


def get_files(paths):
    """Return a list with all the files in paths."""
    sources = []

    for path in paths:
        tmp = [y for x in os.walk(path) for y in glob(
            os.path.join(x[0], '*')) if not os.path.isdir(y)]
        sources.append(tmp)

    return sources


extract_boost()
setup(
    name='pyvinecopulib',
    version='0.0.2',
    install_requires=get_requirements(),
    author='Thomas Nagler and Thibault Vatter',
    author_email='info@vinecopulib.org',
    description='A python interface to vinecopulib',
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/pyvinecopulib/",
    ext_modules=[
        Extension(
            'pyvinecopulib',
            ['src/main.cpp'],
            include_dirs=get_include_paths(),
            depends=get_files(get_include_paths()),
            language='c++'
        ),
    ],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Financial and Insurance Industry',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Other Audience',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Telecommunications Industry',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Programming Language :: C++',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: MIT License'
    ],
    keywords='copula, vines copulas, pair-copulas constructions',
)
