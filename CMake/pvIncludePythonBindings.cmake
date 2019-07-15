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

if(NOT PYVINECOPULIB_PYTHON_OUTPUT_DIRECTORY OR PYVINECOPULIB_PYTHON_OUTPUT_DIRECTORY STREQUAL "")
  set(PYVINECOPULIB_PYTHON_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

if(BUILD_Python_Boost OR BUILD_Python_PyBind)
  include_directories(${PYTHON_INCLUDE_DIRS})
  if(WIN32)
    link_libraries(${PYTHON_LIBRARIES})
  endif()
endif()
