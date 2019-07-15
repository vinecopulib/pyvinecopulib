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
# FLANN - external project needed by PCL.
#-----------------------------------------------------------------------------
if(NOT BUILD_FLANN)
  return()
endif()

# Sanity checks
if(DEFINED FLANN_DIR AND NOT EXISTS ${FLANN_DIR})
  message(FATAL_ERROR "FLANN_DIR variable is defined but corresponds to non-existing directory \"${FLANN_DIR}\".")
endif()

set(version "1.9.1")
set(location "${NIFTK_EP_TARBALL_LOCATION}/flann-${version}.tar.gz")
mpMacroDefineExternalProjectVariables(FLANN ${version} ${location})
set(proj_DEPENDENCIES )

if(NOT DEFINED FLANN_DIR)

  ExternalProject_Add(${proj}
    LIST_SEPARATOR ^^
    PREFIX ${proj_CONFIG}
    SOURCE_DIR ${proj_SOURCE}
    BINARY_DIR ${proj_BUILD}
    INSTALL_DIR ${proj_INSTALL}
    URL ${proj_LOCATION}
    URL_MD5 ${proj_CHECKSUM}
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${EP_COMMON_ARGS}
      -DCMAKE_PREFIX_PATH:PATH=${PYVINECOPULIB_PREFIX_PATH}
      -DUSE_OPENMP:BOOL=${OpenMP_FOUND}
      -DBUILD_CUDA_LIB:BOOL=${PYVINECOPULIB_USE_CUDA}
      -DCUDA_TOOLKIT_ROOT_DIR:PATH=${CUDA_TOOLKIT_ROOT_DIR}
      -DCUDA_ARCH_BIN:STRING=${PYVINECOPULIB_CUDA_ARCH_BIN}
      #-DUSE_MPI:BOOL=${PYVINECOPULIB_USE_MPI} Turning off for now, as you need a parallel version of HDF5.
      -DBUILD_MATLAB_BINDINGS:BOOL=OFF
      -DBUILD_PYTHON_BINDINGS:BOOL=OFF
      -DBUILD_C_BINDINGS:BOOL=OFF
      -DBUILD_TESTS:BOOL=OFF
      -DCMAKE_DEBUG_POSTFIX:STRING=
    CMAKE_CACHE_ARGS
      ${EP_COMMON_CACHE_ARGS}
    CMAKE_CACHE_DEFAULT_ARGS
      ${EP_COMMON_CACHE_DEFAULT_ARGS}
    DEPENDS ${proj_DEPENDENCIES}
  )

  set(FLANN_DIR ${proj_INSTALL})
  set(FLANN_ROOT ${FLANN_DIR})

  set(PYVINECOPULIB_PREFIX_PATH ${proj_INSTALL}^^${PYVINECOPULIB_PREFIX_PATH})
  mitkFunctionInstallExternalCMakeProject(${proj})

  message("SuperBuild loading FLANN from ${FLANN_DIR}")

else()

  mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

endif()
