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

######################################################################
# Configure Dart testing support.  This should be done before any
# message(FATAL_ERROR ...) commands are invoked.
######################################################################

include(${CMAKE_ROOT}/Modules/Dart.cmake)
mark_as_advanced(TCL_TCLSH DART_ROOT)

enable_testing()

if(BUILD_TESTING)
  set(BUILDNAME "PYVINECOPULIB" CACHE STRING "Name of build on the dashboard")
  mark_as_advanced(BUILDNAME)
  configure_file(CMake/CTestCustom.cmake.in ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)
  configure_file(CMake/CTestContinuous.cmake.in ${CMAKE_BINARY_DIR}/CTestContinuous.cmake @ONLY)
endif(BUILD_TESTING)
