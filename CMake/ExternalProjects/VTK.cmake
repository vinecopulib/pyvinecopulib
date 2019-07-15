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
# VTK
#-----------------------------------------------------------------------------
if(NOT BUILD_VTK)
  return()
endif()

# Sanity checks
if(DEFINED VTK_DIR AND NOT EXISTS ${VTK_DIR})
  message(FATAL_ERROR "VTK_DIR variable is defined but corresponds to non-existing directory")
endif()

set(make_webkit_optional)
if( "${VTK_VERSION}" STREQUAL "v6.1.0") # As the patch is only valid for this version.
  set(make_webkit_optional COMMAND ${PATCH_COMMAND} -N -p1 -i ${CMAKE_CURRENT_LIST_DIR}/VTK.patch)
endif()

set(location "https://gitlab.kitware.com/vtk/vtk.git")
mpMacroDefineExternalProjectVariables(VTK ${VTK_VERSION} ${location})
set(proj_DEPENDENCIES )

if(WIN32)
  option(VTK_USE_SYSTEM_FREETYPE OFF)
else(WIN32)
  option(VTK_USE_SYSTEM_FREETYPE ON)
endif(WIN32)
mark_as_advanced(VTK_USE_SYSTEM_FREETYPE)

if(NOT DEFINED VTK_DIR)

  set(additional_cmake_args )

  if(MINGW)
    set(additional_cmake_args
        -DCMAKE_USE_WIN32_THREADS:BOOL=ON
        -DCMAKE_USE_PTHREADS:BOOL=OFF
        -DVTK_USE_VIDEO4WINDOWS:BOOL=OFF # no header files provided by MinGW
        )
  endif()

  ##############################################################################
  # Module selection logic.
  # When running on travis/appveyor, and building a small library that will
  # have a python interface, you will probably want the smallest build possible.
  # This will promote the idea that a small python extension should be as small
  # as possible, and just provide a few re-usable algorithms.
  #
  # So, as an example, the next if block shows how to turn most default
  # modules OFF, and then just turn on one set of algorithms.
  ##############################################################################
  if(BUILD_Python_Boost OR BUILD_Python_PyBind)

    # So, a minimum build, plus one module.
    list(APPEND additional_cmake_args
      -DVTK_Group_Rendering:BOOL=OFF
      -DVTK_Group_StandAlone:BOOL=OFF
      -DModule_vtkFiltersGeneral:BOOL=ON
    )
  else()

    # Otherwise, we will build all default modules.

    message("Building mostly default VTK modules")

    # And then this additionally turns on vtkRenderingExternal
    # which is used for QML rendering, and valid from v7.1.0 onwards.
    # Look in pvAddVTK, as we either build v6.1.0, or v8.2.0.
    if( NOT "${VTK_VERSION}" STREQUAL "v6.1.0")
      list(APPEND additional_cmake_args
        -DModule_vtkRenderingExternal:BOOL=ON
      )
    endif()

  endif()

  # Optionally enable memory leak checks for any objects derived from vtkObject. This
  # will force unit tests to fail if they have any of these memory leaks.
  option(PYVINECOPULIB_VTK_DEBUG_LEAKS OFF)
  mark_as_advanced(PYVINECOPULIB_VTK_DEBUG_LEAKS)
  list(APPEND additional_cmake_args
       -DVTK_DEBUG_LEAKS:BOOL=${MITK_VTK_DEBUG_LEAKS}
      )

  list(APPEND additional_cmake_args
       -DVTK_WRAP_PYTHON:BOOL=OFF
       -DVTK_WINDOWS_PYTHON_DEBUGGABLE:BOOL=OFF
      )

  if(Qt5_DIR AND PYVINECOPULIB_USE_QT)
    list(APPEND additional_cmake_args
        -DVTK_QT_VERSION:STRING=5
        -DVTK_Group_Qt:BOOL=ON
        -DVTK_INSTALL_NO_QT_PLUGIN:BOOL=ON
     )
  endif()

  if(CTEST_USE_LAUNCHERS)
    list(APPEND additional_cmake_args
      "-DCMAKE_PROJECT_${proj}_INCLUDE:FILEPATH=${CMAKE_ROOT}/Modules/CTestUseLaunchers.cmake"
    )
  endif()

  if(APPLE)
    list(APPEND additional_cmake_args
        -DVTK_REQUIRED_OBJCXX_FLAGS:STRING=""
        )
  endif(APPLE)

  ExternalProject_Add(${proj}
    LIST_SEPARATOR ^^
    PREFIX ${proj_CONFIG}
    SOURCE_DIR ${proj_SOURCE}
    BINARY_DIR ${proj_BUILD}
    INSTALL_DIR ${proj_INSTALL}
    GIT_REPOSITORY ${proj_LOCATION}
    GIT_TAG ${proj_VERSION}
    UPDATE_COMMAND ${GIT_EXECUTABLE} checkout ${proj_VERSION}
    PATCH_COMMAND ${make_webkit_optional}
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
        ${EP_COMMON_ARGS}
        -DCMAKE_PREFIX_PATH:PATH=${PYVINECOPULIB_PREFIX_PATH}
        -DVTK_WRAP_TCL:BOOL=OFF
        -DVTK_WRAP_PYTHON:BOOL=OFF
        -DVTK_WRAP_JAVA:BOOL=OFF
        -DVTK_USE_SYSTEM_FREETYPE:BOOL=${VTK_USE_SYSTEM_FREETYPE}
        -DVTK_LEGACY_REMOVE:BOOL=ON
        -DVTK_MAKE_INSTANTIATORS:BOOL=ON
        -DVTK_USE_CXX11_FEATURES:BOOL=ON
        -DVTK_RENDERING_BACKEND:STRING=${VTK_BACKEND}
        ${additional_cmake_args}
    CMAKE_CACHE_ARGS
      ${EP_COMMON_CACHE_ARGS}
    CMAKE_CACHE_DEFAULT_ARGS
      ${EP_COMMON_CACHE_DEFAULT_ARGS}
    DEPENDS ${proj_DEPENDENCIES}
  )

  set(VTK_DIR ${proj_INSTALL})
  set(PYVINECOPULIB_PREFIX_PATH ${proj_INSTALL}^^${PYVINECOPULIB_PREFIX_PATH})
  mitkFunctionInstallExternalCMakeProject(${proj})

  message("SuperBuild loading VTK from ${VTK_DIR}")

else()

  mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

endif()
