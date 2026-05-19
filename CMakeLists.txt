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

# Replace upstream's `-march=native` with v3 for distributed wheels so they run
# on every CPU since ~2013. Editable builds keep `-march=native`.
if(DEFINED ENV{CIBUILDWHEEL})
  string(REGEX REPLACE "(^| )-march=native( |$)" " -march=x86-64-v3 "
                       CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

set(vinecopulib_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/vinecopulib/include)
set(wdm_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/wdm/include)
set(kde1d_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/kde1d/include)
set(pyvinecopulib_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/src/include)

set(BUILD_TESTING OFF)

if(DEFINED ENV{EIGEN3_INCLUDE_DIR})
  set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
endif()

if(DEFINED ENV{Boost_INCLUDE_DIR})
  set(Boost_INCLUDE_DIRS $ENV{Boost_INCLUDE_DIR})
endif()

include(lib/vinecopulib/cmake/findDependencies.cmake REQUIRED)

find_package(
  Python 3.8 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
message(STATUS "Python_EXECUTABLE= ${Python_EXECUTABLE}")

find_package(nanobind CONFIG REQUIRED)

include(lib/vinecopulib/cmake/printInfo.cmake REQUIRED)

include_directories(SYSTEM ${external_includes})
include_directories(${vinecopulib_INCLUDE_DIRS} ${kde1d_INCLUDE_DIRS}
                    ${pyvinecopulib_INCLUDE_DIRS})

# Regenerate src/include/docstr.hpp from the C++ headers via libclang.
# Gitignored: the build is the single source of truth.
execute_process(
  COMMAND "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/find_libclang.py"
          --search-paths "${CMAKE_PREFIX_PATH}"
  OUTPUT_VARIABLE LIBCLANG_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT LIBCLANG_PATH OR NOT EXISTS "${LIBCLANG_PATH}")
  message(
    FATAL_ERROR
      "Could not locate libclang. Install via 'pip install libclang' "
      "(already in [build-system] requires) or include python-clang in "
      "your conda env.")
endif()
message(STATUS "Found libclang: ${LIBCLANG_PATH}")

# Both .hpp and .ipp are passed so the script can extract symbols defined in
# .ipp files; the glue TU only #includes .hpp to avoid Windows redefinitions
# (.ipp lacks include guards). Trailing `<key>.<ext>` (no wildcard) prevents
# matching unrelated files like `tools_stats_ghalton.hpp`.
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
    "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/generate_docstring.py"
    -std "c++${CMAKE_CXX_STANDARD}" -I
    "${CMAKE_SOURCE_DIR}/lib/vinecopulib/include" -I
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

nanobind_add_module(pyvinecopulib_ext STABLE_ABI LTO NOMINSIZE NB_STATIC
                    src/pyvinecopulib_ext.cpp)

add_dependencies(pyvinecopulib_ext generate_docstr)

target_compile_definitions(pyvinecopulib_ext
                           PRIVATE VERSION_INFO=\"${PROJECT_VERSION}\")

# Generate per-subpackage __init__.pyi stubs. The sklearn stub is OPTIONAL (only
# produced when scikit-learn is available at build time).
add_custom_command(
  TARGET pyvinecopulib_ext
  POST_BUILD
  COMMAND
    "${Python_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/generate_stubs.py"
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

install(TARGETS pyvinecopulib_ext LIBRARY DESTINATION pyvinecopulib)

# Install stubs + py.typed explicitly (gitignored, so git-ls-files skips them).
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
