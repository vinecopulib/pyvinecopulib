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

if(NOT BUILDING_GUIS)

  if (NOT PYVINECOPULIB_PYTHON_MODULE_NAME OR PYVINECOPULIB_PYTHON_MODULE_NAME STREQUAL "")
    set(PYVINECOPULIB_PYTHON_MODULE_NAME "" CACHE STRING "Basename for python module." FORCE)
  endif()

  set(PYVINECOPULIB_PYTHON_SRCS
    pvLibPython.cpp
  )
  add_subdirectory(pybind11)
  link_libraries(${ALL_LIBRARIES})
  pybind11_add_module(${PYVINECOPULIB_PYTHON_MODULE_NAME} ${PYVINECOPULIB_PYTHON_SRCS})
  if(APPLE)
    set_target_properties(${PYVINECOPULIB_PYTHON_MODULE_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif()
  install(TARGETS ${PYVINECOPULIB_PYTHON_MODULE_NAME} DESTINATION ${PYVINECOPULIB_INSTALL_LIB_DIR})

endif()
