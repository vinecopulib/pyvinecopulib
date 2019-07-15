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

if(MSVC)
  set(variables
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
  )
  if(BUILD_SHARED_LIBS)
    message(STATUS "Forcing MSVC to use dynamic runtime variables.")
    foreach(variable ${variables})
      if(${variable} MATCHES "/MT")
        string(REGEX REPLACE "/MT" "/MD" ${variable} "${${variable}}")
      endif()
    endforeach()
  else()
    message(STATUS "Forcing MSVC to use static runtime variables.")
    foreach(variable ${variables})
      if(${variable} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${variable} "${${variable}}")
      endif()
    endforeach()
  endif()
  message(STATUS "Initial build flags:")
  foreach(variable ${variables})
    message(STATUS "  '${variable}': ${${variable}}")
  endforeach()
  message(STATUS "")
endif()
