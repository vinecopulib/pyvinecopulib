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
# ArrayFire
#-----------------------------------------------------------------------------
if(NOT BUILD_ArrayFire)
  return()
endif()

# Sanity checks
if(DEFINED ArrayFire_DIR AND NOT EXISTS ${ArrayFire_DIR})
  message(FATAL_ERROR "ArrayFire_DIR variable is defined but corresponds to non-existing directory \"${ArrayFire_DIR}\".")
endif()

set(ArrayFire_VERSION "dc38ef1329")
set(location "https://github.com/arrayfire/arrayfire.git")
mpMacroDefineExternalProjectVariables(ArrayFire ${ArrayFire_VERSION} ${location})
set(proj_DEPENDENCIES )

if(NOT DEFINED ArrayFire_DIR)

  ExternalProject_Add(${proj}
    LIST_SEPARATOR ^^
    PREFIX ${proj_CONFIG}
    SOURCE_DIR ${proj_SOURCE}
    BINARY_DIR ${proj_BUILD}
    INSTALL_DIR ${proj_INSTALL}
    GIT_REPOSITORY ${proj_LOCATION}
    GIT_TAG ${proj_VERSION}
    GIT_SHALLOW ON
    UPDATE_COMMAND ${GIT_EXECUTABLE} checkout ${proj_VERSION}
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${EP_COMMON_ARGS}
      -DCMAKE_PREFIX_PATH:PATH=${PYVINECOPULIB_PREFIX_PATH}
      -DAF_BUILD_CPU:BOOL=ON
      -DAF_BUILD_CUDA:BOOL=OFF
      -DAF_BUILD_OPENCL:BOOL=OFF
    CMAKE_CACHE_ARGS
      ${EP_COMMON_CACHE_ARGS}
    CMAKE_CACHE_DEFAULT_ARGS
      ${EP_COMMON_CACHE_DEFAULT_ARGS}
    DEPENDS ${proj_DEPENDENCIES}
  )

  set(ArrayFire_DIR ${proj_INSTALL})
  set(PYVINECOPULIB_PREFIX_PATH ${proj_INSTALL}^^${PYVINECOPULIB_PREFIX_PATH})
  mitkFunctionInstallExternalCMakeProject(${proj})

  message("SuperBuild loading ArrayFire from ${ArrayFire_DIR}")

else()

  mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

endif()
