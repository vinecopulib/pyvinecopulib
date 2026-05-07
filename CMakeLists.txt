cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(
  pyvinecopulib
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

if(NOT SKBUILD)
  message(
    WARNING
      "\
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

include(lib/vinecopulib/cmake/options.cmake REQUIRED)

include(lib/vinecopulib/cmake/compilerDefOpt.cmake REQUIRED)

# wdm_INCLUDE_DIRS is lib/wdm/include
set(vinecopulib_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/vinecopulib/include)
set(wdm_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/wdm/include)
set(kde1d_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/kde1d/include)
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

include(lib/vinecopulib/cmake/findDependencies.cmake REQUIRED)

# Try to import all Python components potentially needed by nanobind
find_package(
  Python 3.8 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
message(STATUS "Python_EXECUTABLE= ${Python_EXECUTABLE}")

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

include(lib/vinecopulib/cmake/printInfo.cmake REQUIRED)

include_directories(SYSTEM ${external_includes})
include_directories(${vinecopulib_INCLUDE_DIRS} ${kde1d_INCLUDE_DIRS}
                    ${pyvinecopulib_INCLUDE_DIRS})

# ==============================================================================
# Generated docstrings: src/include/docstr.hpp
#
# Regenerated on every configure/build via libclang parsing of the C++ headers
# in lib/. The output is gitignored — the build is the single source of truth,
# eliminating the need for a regenerate-and-diff drift gate. libclang is
# declared as a build-system requirement in pyproject.toml, so isolated builds
# pull a known-good version automatically.
# ==============================================================================

execute_process(
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/find_libclang.py"
  OUTPUT_VARIABLE LIBCLANG_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT LIBCLANG_PATH OR NOT EXISTS "${LIBCLANG_PATH}")
  message(
    FATAL_ERROR
      "Could not locate libclang. Required at build time to regenerate "
      "src/include/docstr.hpp. Install via 'pip install libclang' (already "
      "in [build-system] requires for isolated builds) or include "
      "python-clang in your conda env.")
endif()
message(STATUS "Found libclang: ${LIBCLANG_PATH}")

# Headers parsed for docstring extraction. Mirrors the substring filter used by
# the previous scripts/generate_metadata.py (class/controls/family/structure
# /stats/wdm/kde1d × .hpp/.ipp). Globs run at CMake configure time.
file(
  GLOB_RECURSE
  DOCSTR_INPUT_HEADERS
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*class*.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*class*.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*controls*.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*controls*.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*family*.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*family*.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*structure*.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*structure*.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*stats*.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*stats*.ipp"
  "${CMAKE_SOURCE_DIR}/lib/wdm/include/*wdm*.hpp"
  "${CMAKE_SOURCE_DIR}/lib/kde1d/include/*kde1d*.hpp"
  "${CMAKE_SOURCE_DIR}/lib/kde1d/include/*stats*.hpp")

set(DOCSTR_OUTPUT "${CMAKE_SOURCE_DIR}/src/include/docstr.hpp")

add_custom_command(
  OUTPUT "${DOCSTR_OUTPUT}"
  COMMAND
    "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/generate_docstring.py"
    -I "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include" -I
    "${CMAKE_SOURCE_DIR}/lib/wdm/include" -I
    "${CMAKE_SOURCE_DIR}/lib/kde1d/include" -isystem "${EIGEN3_INCLUDE_DIR}"
    -isystem "${Boost_INCLUDE_DIRS}" -output "${DOCSTR_OUTPUT}" -library_file
    "${LIBCLANG_PATH}" ${DOCSTR_INPUT_HEADERS}
  DEPENDS ${DOCSTR_INPUT_HEADERS}
          "${CMAKE_SOURCE_DIR}/scripts/generate_docstring.py"
          "${CMAKE_SOURCE_DIR}/scripts/find_libclang.py"
  COMMENT "Generating docstr.hpp from C++ headers via libclang"
  VERBATIM)

add_custom_target(generate_docstr DEPENDS "${DOCSTR_OUTPUT}")

nanobind_add_module(
  pyvinecopulib_ext
  # Target the stable ABI for Python 3.12+, which reduces the number of binary
  # wheels that must be built. This does nothing on older Python versions
  STABLE_ABI
  # Enable link-time optimization for the extension
  LTO
  # Don’t perform optimizations to minimize binary size.
  NOMINSIZE
  # Build libnanobind statically and merge it into the extension (which itself
  # remains a shared library)
  #
  # If your project builds multiple extensions, you can replace this flag by
  # NB_SHARED to conserve space by reusing a shared libnanobind across libraries
  NB_STATIC
  src/pyvinecopulib_ext.cpp)

# Ensure docstr.hpp exists before pyvinecopulib_ext.cpp is compiled.
add_dependencies(pyvinecopulib_ext generate_docstr)

target_compile_definitions(pyvinecopulib_ext
                           PRIVATE VERSION_INFO=\"${PROJECT_VERSION}\")

# ==============================================================================
# Generated stubs: src/pyvinecopulib/__init__.pyi
#
# After the extension is built, stage the package + the just-built .so in a
# tempdir and run scripts/generate_stubs.py against it. Output is gitignored.
# ==============================================================================

add_custom_command(
  TARGET pyvinecopulib_ext
  POST_BUILD
  COMMAND
    "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/build_stubs.py"
    --source-pkg-dir "${CMAKE_SOURCE_DIR}/src/pyvinecopulib" --ext-so
    "$<TARGET_FILE:pyvinecopulib_ext>" --output
    "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/__init__.pyi" --py-typed
    "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/py.typed"
  COMMENT "Generating __init__.pyi from the just-built extension"
  VERBATIM)

# Install directive for scikit-build-core
install(TARGETS pyvinecopulib_ext LIBRARY DESTINATION pyvinecopulib)
