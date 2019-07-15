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

#-----------------------------------------------------------------------------
# gflags
#-----------------------------------------------------------------------------
if(NOT BUILD_gflags)
  return()
endif()

# Sanity checks
if(DEFINED gflags_DIR AND NOT EXISTS ${gflags_DIR})
  message(FATAL_ERROR "gflags_DIR variable is defined but corresponds to non-existing directory \"${gflags_ROOT}\".")
endif()

set(gflags_VERSION "46f73f88b1")
set(location "https://github.com/gflags/gflags.git")
mpMacroDefineExternalProjectVariables(gflags ${gflags_VERSION} ${location})
set(proj_DEPENDENCIES )

if(NOT DEFINED gflags_DIR)

  ExternalProject_Add(${proj}
    LIST_SEPARATOR ^^
    PREFIX ${proj_CONFIG}
    SOURCE_DIR ${proj_SOURCE}
    BINARY_DIR ${proj_BUILD}
    INSTALL_DIR ${proj_INSTALL}
    GIT_REPOSITORY ${proj_LOCATION}
    GIT_TAG ${proj_VERSION}
    UPDATE_COMMAND ${GIT_EXECUTABLE} checkout ${proj_VERSION}
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${EP_COMMON_ARGS}
      -DCMAKE_PREFIX_PATH:PATH=${PYVINECOPULIB_PREFIX_PATH}
    CMAKE_CACHE_ARGS
      ${EP_COMMON_CACHE_ARGS}
    CMAKE_CACHE_DEFAULT_ARGS
      ${EP_COMMON_CACHE_DEFAULT_ARGS}
    DEPENDS ${proj_DEPENDENCIES}
  )

  set(gflags_DIR ${proj_INSTALL})
  set(gflags_LIBRARY_DIR ${gflags_DIR}/lib)     # gflags-config.cmake doesn't export this.
  set(gflags_INCLUDE_DIR ${gflags_DIR}/include)

  mitkFunctionInstallExternalCMakeProject(${proj})

  message("SuperBuild loading gflags from ${gflags_DIR}.")

else(NOT DEFINED gflags_DIR)

  mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

endif(NOT DEFINED gflags_DIR)
