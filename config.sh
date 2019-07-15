#/*============================================================================
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
#============================================================================*/

function pre_build {
  echo "Starting pre_build."

  # Debug info.
  pwd
  echo "PATH=$PATH"
  python --version
  cmake --version

  if [ -n "$BUILD_DEPENDS" ]; then
    pip install $BUILD_DEPENDS
  fi

  if [ -n "$IS_OSX" ]; then
    echo "pre_build is on Mac."
  else
    echo "pre_build is on Linux."

    #############################################################################################
    # IMPORTANT: Look in .travis.yml. Decide if your project requires DO_PYTHON_BUILD to be true.
    #############################################################################################
    #
    # if DO_PYTHON_BUILD = true
    #
    #   i.e. you want to build C++ and then Python wheels. Therefore, this function is running inside
    #   a manylinux docker image. Any library installations here should use Centos based commands
    #   like 'sudo yum install ...'.
    #   OR
    #   You should build your own docker image, specify that URL in .travis.yml,
    #   and then you won't need any 'sudo yum install ...' type commands here.
    #
    # if DO_PYTHON_BUILD = false
    #
    #   then this function is running inside the main travis VM, most likely Ubuntu,
    #   so you should set up dependencies in top-level .travis.yml, or try some
    #   'sudo apt-get install' type commands here. But its more obvious if you
    #   put them in the top level .travis.yml file.

  fi

  # Run the actual C++ build.
  source travis_cmake_build.sh
  cmake_build

  echo "Finished pre_build."
}

function build_wheel {
  # Don't remove this function, or switch back to pip, as pip
  # doesn't get on well with versioneer.py
  # https://github.com/warner/python-versioneer/issues/121
  build_bdist_wheel $@
}

function run_tests {
  echo "Starting run_tests."
  pwd
  cd ..
  # Uncomment this when you have some tests.
  # If no tests are found, this fails with non-zero error code.
  # python -m pytest  -v -s Testing/
  echo "Finished run_tests."
}
