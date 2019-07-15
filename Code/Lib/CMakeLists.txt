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

set(PYVINECOPULIB_LIBRARY_NAME pyvinecopulib)

set(PYVINECOPULIB_LIBRARY_SRCS
  pvBasicTypes.cpp
  pvMyFunctions.cpp
  pvException.cpp
)

set(PYVINECOPULIB_LIBRARY_HDRS
  pvExceptionMacro.h
)

add_library(${PYVINECOPULIB_LIBRARY_NAME} ${PYVINECOPULIB_LIBRARY_HDRS} ${PYVINECOPULIB_LIBRARY_SRCS})

if(BUILD_SHARED_LIBS)
  if(WIN32)
    add_definitions(-DPYVINECOPULIB_WINDOWS_EXPORT)
  endif(WIN32)
endif(BUILD_SHARED_LIBS)

target_link_libraries(${PYVINECOPULIB_LIBRARY_NAME} PRIVATE ${ALL_THIRD_PARTY_LIBRARIES})

PYVINECOPULIB_INSTALL_HEADERS()
PYVINECOPULIB_INSTALL_LIBRARY(${PYVINECOPULIB_LIBRARY_NAME})

