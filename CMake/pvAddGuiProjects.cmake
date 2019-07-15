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

# Todo: Make macro to simplify this config.

option(BUILD_QtVTKDemo "Build QtVTKDemo Gui." OFF)
option(BUILD_QMLVTKDemo "Build QMLVTKDemo Gui." OFF)
option(BUILD_QOpenGLDemo "Build QOpenGLDemo Gui." OFF)

if(BUILD_PYTHON_BINDINGS AND BUILD_QtVTKDemo)
  set(BUILD_QtVTKDemo OFF CACHE BOOL "Build QtVTKDemo Gui." FORCE)
  message("Forcing BUILD_QtVTKDemo to OFF as you want a python module.")
endif()

if(BUILD_PYTHON_BINDINGS AND BUILD_QMLVTKDemo)
  set(BUILD_QMLVTKDemo OFF CACHE BOOL "Build QMLVTKDemo Gui." FORCE)
  message("Forcing BUILD_QMLVTKDemo to OFF as you want a python module.")
endif()

if(BUILD_PYTHON_BINDINGS AND BUILD_QOpenGLDemo)
  set(BUILD_QOpenGLDemo OFF CACHE BOOL "Build QOpenGLDemo Gui." FORCE)
  message("Forcing BUILD_QOpenGLDemo to OFF as you want a python module.")
endif()

option(PYVINECOPULIB_USE_QT "Use Qt." OFF)
mark_as_advanced(PYVINECOPULIB_USE_QT) # Qt gets baked into VTK, so really developers should not fiddle with this.

if(BUILD_QtVTKDemo AND NOT PYVINECOPULIB_USE_QT)
  set(PYVINECOPULIB_USE_QT ON CACHE BOOL "Use Qt." FORCE)
  message("Forcing PYVINECOPULIB_USE_QT to ON due to BUILD_QtVTKDemo being ON.")
endif()

if(BUILD_QtVTKDemo AND NOT BUILD_VTK)
  set(BUILD_VTK ON CACHE BOOL "Build VTK." FORCE)
  message("Forcing BUILD_VTK to ON due to BUILD_QtVTKDemo being ON.")
endif()

if(BUILD_QMLVTKDemo AND NOT PYVINECOPULIB_USE_QT)
  set(PYVINECOPULIB_USE_QT ON CACHE BOOL "Use Qt." FORCE)
  message("Forcing PYVINECOPULIB_USE_QT to ON due to BUILD_QMLVTKDemo being ON.")
endif()

if(BUILD_QMLVTKDemo AND NOT BUILD_VTK)
  set(BUILD_VTK ON CACHE BOOL "Build VTK." FORCE)
  message("Forcing BUILD_VTK to ON due to BUILD_QMLVTKDemo being ON.")
endif()

if(BUILD_QMLVTKDemo AND NOT BUILD_SHARED_LIBS)
  message("Forcing BUILD_SHARED_LIBS to ON due to BUILD_QMLVTKDemo being ON.")
  set(BUILD_SHARED_LIBS ON CACHE BOOL "Build Shared Libraries" FORCE)
endif()

if(BUILD_QOpenGLDemo AND NOT PYVINECOPULIB_USE_QT)
  set(PYVINECOPULIB_USE_QT ON CACHE BOOL "Use Qt." FORCE)
  message("Forcing PYVINECOPULIB_USE_QT to ON due to BUILD_QOpenGLDemo being ON.")
endif()

######################################################################
# This is a list of all known GUI apps. Used for things like
# creating Mac OSX bundles, and each command line app gets copied
# into each OSX bundle.
######################################################################
if (BUILD_QtVTKDemo)
  list(APPEND _known_apps QtVTKDemo)
  set(BUILDING_GUIS ON)
endif()
if (BUILD_QMLVTKDemo)
  list(APPEND _known_apps QMLVTKDemo)
  set(BUILDING_GUIS ON)
endif()
if (BUILD_QOpenGLDemo)
  list(APPEND _known_apps QOpenGLDemo)
  set(BUILDING_GUIS ON)
endif()
set_property(GLOBAL PROPERTY PYVINECOPULIB_GUI_APPS ${_known_apps})
