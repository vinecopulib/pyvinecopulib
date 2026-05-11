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
  COMMAND
    "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/build/find_libclang.py"
    --search-paths "${CMAKE_PREFIX_PATH}"
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

# Headers parsed for docstring extraction. Both .hpp and .ipp are passed: the
# script tracks both in its include-file map (so symbols defined in .ipp
# implementation files are extracted), but only #includes .hpp files in the glue
# translation unit (.ipp content reaches libclang transitively via each .hpp's
# bottom include, and double inclusion would break Windows since .ipp has no
# include guards). Patterns end with `<key>.<ext>` (no trailing wildcard) to
# avoid matching unrelated files like `tools_stats_ghalton.hpp`.
file(
  GLOB_RECURSE
  DOCSTR_INPUT_HEADERS
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*class.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*class.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*controls.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*controls.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*family.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*family.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*structure.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*structure.ipp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*stats.hpp"
  "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include/*stats.ipp"
  "${CMAKE_SOURCE_DIR}/lib/wdm/include/*wdm.hpp"
  "${CMAKE_SOURCE_DIR}/lib/kde1d/include/*kde1d.hpp"
  "${CMAKE_SOURCE_DIR}/lib/kde1d/include/*stats.hpp")

set(DOCSTR_OUTPUT "${CMAKE_SOURCE_DIR}/src/include/docstr.hpp")

add_custom_command(
  OUTPUT "${DOCSTR_OUTPUT}"
  COMMAND
    "${Python_EXECUTABLE}"
    "${CMAKE_SOURCE_DIR}/scripts/build/generate_docstring.py" -std
    "c++${CMAKE_CXX_STANDARD}" -I "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include"
    -I "${CMAKE_SOURCE_DIR}/lib/wdm/include" -I
    "${CMAKE_SOURCE_DIR}/lib/kde1d/include" -isystem "${EIGEN3_INCLUDE_DIR}"
    -isystem "${Boost_INCLUDE_DIRS}" -output "${DOCSTR_OUTPUT}" -library_file
    "${LIBCLANG_PATH}" ${DOCSTR_INPUT_HEADERS}
  DEPENDS ${DOCSTR_INPUT_HEADERS}
          "${CMAKE_SOURCE_DIR}/scripts/build/generate_docstring.py"
          "${CMAKE_SOURCE_DIR}/scripts/build/find_libclang.py"
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
# Generated stubs: one __init__.pyi per public (sub)package.
#
# After the extension is built, scripts/build/generate_stubs.py stages the
# package + the just-built .so in a tempdir and writes a .pyi for each requested
# module. Outputs are gitignored.
#
# pyvinecopulib.sklearn is optional (requires scikit-learn on the build
# machine); the generator emits a warning and skips it if the import fails.
# ==============================================================================

add_custom_command(
  TARGET pyvinecopulib_ext
  POST_BUILD
  COMMAND
    "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/build/generate_stubs.py"
    --source-pkg-dir "${CMAKE_SOURCE_DIR}/src/pyvinecopulib" --ext-so
    "$<TARGET_FILE:pyvinecopulib_ext>" --py-typed
    "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/py.typed" --module-output
    "pyvinecopulib:${CMAKE_SOURCE_DIR}/src/pyvinecopulib/__init__.pyi"
    --module-output
    "pyvinecopulib.core:${CMAKE_SOURCE_DIR}/src/pyvinecopulib/core/__init__.pyi"
    --module-output
    "pyvinecopulib.families:${CMAKE_SOURCE_DIR}/src/pyvinecopulib/families/__init__.pyi"
    --module-output
    "pyvinecopulib.utils:${CMAKE_SOURCE_DIR}/src/pyvinecopulib/utils/__init__.pyi"
    --module-output
    "pyvinecopulib.sklearn:${CMAKE_SOURCE_DIR}/src/pyvinecopulib/sklearn/__init__.pyi"
  COMMENT
    "Generating per-subpackage __init__.pyi stubs from the just-built extension"
  VERBATIM)

# Install directive for scikit-build-core
install(TARGETS pyvinecopulib_ext LIBRARY DESTINATION pyvinecopulib)

# Explicitly install the generated stubs and py.typed marker into the wheel. All
# are gitignored build artifacts; scikit-build-core's git-ls-files-based package
# discovery would otherwise skip them. The sklearn stub is OPTIONAL: it's only
# produced if scikit-learn is available on the build machine, so the install
# rule must tolerate its absence.
install(FILES "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/__init__.pyi"
              "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/py.typed"
        DESTINATION pyvinecopulib)
install(FILES "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/core/__init__.pyi"
        DESTINATION pyvinecopulib/core)
install(FILES "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/families/__init__.pyi"
        DESTINATION pyvinecopulib/families)
install(FILES "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/utils/__init__.pyi"
        DESTINATION pyvinecopulib/utils)
install(
  FILES "${CMAKE_SOURCE_DIR}/src/pyvinecopulib/sklearn/__init__.pyi"
  DESTINATION pyvinecopulib/sklearn
  OPTIONAL)
