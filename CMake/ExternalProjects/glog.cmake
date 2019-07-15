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
# glog
#-----------------------------------------------------------------------------
if(NOT BUILD_glog)
  return()
endif()

# Sanity checks
if(DEFINED glog_DIR AND NOT EXISTS ${glog_DIR})
  message(FATAL_ERROR "glog_DIR variable is defined but corresponds to non-existing directory \"${glog_ROOT}\".")
endif()

set(glog_VERSION "a6a166db06")
set(location "https://github.com/google/glog.git")
mpMacroDefineExternalProjectVariables(glog ${glog_VERSION} ${location})
set(proj_DEPENDENCIES gflags)

if(NOT DEFINED glog_DIR)

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
      -Dgflags_DIR:PATH=${gflags_DIR}/lib/cmake/gflags
      -DCMAKE_SHARED_LINKER_FLAGS:STRING="-L${gflags_LIBRARY_DIR}" # gflags-config.cmake doesn't export this.
      -DCMAKE_EXE_LINKER_FLAGS:STRING="-L${gflags_LIBRARY_DIR}"    # gflags-config.cmake doesn't export this.
    CMAKE_CACHE_ARGS
      ${EP_COMMON_CACHE_ARGS}
    CMAKE_CACHE_DEFAULT_ARGS
      ${EP_COMMON_CACHE_DEFAULT_ARGS}
    DEPENDS ${proj_DEPENDENCIES}
  )

  set(glog_DIR ${proj_INSTALL})
  set(glog_INCLUDE_DIR ${glog_DIR}/include)
  set(glog_LIBRARY_DIR ${glog_DIR}/lib)
  mitkFunctionInstallExternalCMakeProject(${proj})

  message("SuperBuild loading glog from ${glog_DIR}.")

else(NOT DEFINED glog_DIR)

  mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

endif(NOT DEFINED glog_DIR)

