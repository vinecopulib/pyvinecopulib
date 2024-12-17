cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)

project(
  pyvinecopulib
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX
)

if (NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

include(lib/vinecopulib/cmake/options.cmake               REQUIRED)

include(lib/vinecopulib/cmake/compilerDefOpt.cmake        REQUIRED)

# wdm_INCLUDE_DIRS is lib/wdm/include
set(vinecopulib_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/vinecopulib/include)
set(wdm_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/wdm/include)
set(pyvinecopulib_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/src/include)

# BUILD_TESTING is turned off
set(BUILD_TESTING OFF)

# Look for Eigen3 and Boost
if(DEFINED ENV{EIGEN3_INCLUDE_DIR})
  set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
endif()

if(DEFINED ENV{Boost_INCLUDE_DIR})
  set(Boost_INCLUDE_DIRS $ENV{Boost_INCLUDE_DIR})
endif()

include(lib/vinecopulib/cmake/findDependencies.cmake      REQUIRED)

# Try to import all Python components potentially needed by nanobind
find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
message(STATUS "Python_EXECUTABLE= ${Python_EXECUTABLE}")

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

include(lib/vinecopulib/cmake/printInfo.cmake             REQUIRED)

include_directories(SYSTEM ${external_includes})
include_directories(${vinecopulib_INCLUDE_DIRS} ${pyvinecopulib_INCLUDE_DIRS})

nanobind_add_module(
  pyvinecopulib_ext
  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  STABLE_ABI

  # Enable link-time optimization for the extension
  LTO

  # Build libnanobind statically and merge it into the
  # extension (which itself remains a shared library)
  #
  # If your project builds multiple extensions, you can
  # replace this flag by NB_SHARED to conserve space by
  # reusing a shared libnanobind across libraries
  NB_STATIC

  src/pyvinecopulib_ext.cpp
)

target_compile_definitions(pyvinecopulib_ext PRIVATE VERSION_INFO=\"${PROJECT_VERSION}\")

# Install directive for scikit-build-core
install(TARGETS pyvinecopulib_ext
    LIBRARY DESTINATION pyvinecopulib
)