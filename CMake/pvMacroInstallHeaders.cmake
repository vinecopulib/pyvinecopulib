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

macro(PYVINECOPULIB_INSTALL_HEADERS)

  set(ARGS ${ARGN})
  list(FIND ARGS DESTINATION _destination_index)
  if(_destination_index GREATER -1)
    message(SEND_ERROR "PYVINECOPULIB_INSTALL_HEADERS macro must not be called with a DESTINATION parameter.")
  else()

    if(NOT BUILDING_GUIS)
      file(GLOB_RECURSE HEADERS *.h)
      install(FILES ${HEADERS} DESTINATION ${PYVINECOPULIB_INSTALL_INC_DIR} COMPONENT HEADERS)
    endif()

  endif()
endmacro()
