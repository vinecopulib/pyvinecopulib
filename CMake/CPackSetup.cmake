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

############################################################################
# First, set the generator variable
############################################################################

if(NOT CPACK_GENERATOR)
  if(WIN32)
  
    find_program(NSIS_MAKENSIS NAMES makensis
      PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\NSIS]
      DOC "Where is makensis.exe located"
      )

    if(NOT NSIS_MAKENSIS)
      set(CPACK_GENERATOR ZIP)
    else()
      set(CPACK_GENERATOR "NSIS")
    endif(NOT NSIS_MAKENSIS)
    
    set(CPACK_SOURCE_GENERATOR ZIP)
    
  else()
  
    if(APPLE)
      set(CPACK_GENERATOR "DragNDrop")
    else()
      set(CPACK_GENERATOR TBZ2)
    endif()
    
    set(CPACK_SOURCE_GENERATOR TBZ2)
    
  endif()
endif(NOT CPACK_GENERATOR)

############################################################################
# This bit came from MITK (http://www.mitk.org). Don't know if we need it. Yes we do.
############################################################################
# Should apply only to windows. Note however, that the debug crt is
# non-redistributable! a debug-package cannot be made available to the public.
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(CMAKE_INSTALL_DEBUG_LIBRARIES ON)
endif()
include(InstallRequiredSystemLibraries)

############################################################################
# This bit came from various CMIC packages - START
############################################################################
if (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
  set (CMAKE_SYSTEM_PROCESSOR "x86")
endif ()
if(NOT DEFINED CPACK_SYSTEM_NAME)
  set(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif()
if(${CPACK_SYSTEM_NAME} MATCHES Windows)
  if(CMAKE_CL_64)
    set(CPACK_SYSTEM_NAME Win64-${CMAKE_SYSTEM_PROCESSOR})
  else()
    set(CPACK_SYSTEM_NAME Win32-${CMAKE_SYSTEM_PROCESSOR})
  endif()
endif()

if(${CPACK_SYSTEM_NAME} MATCHES Darwin AND CMAKE_OSX_ARCHITECTURES)
  list(LENGTH CMAKE_OSX_ARCHITECTURES _length)
  if(_length GREATER 1)
    set(CPACK_SYSTEM_NAME Darwin-Universal)
  else()
    set(CPACK_SYSTEM_NAME Darwin-${CMAKE_OSX_ARCHITECTURES})
  endif()
endif()
############################################################################
# This bit came from various CMIC packages - END
############################################################################

############################################################################
# The main setting of CPack settings that are independent of generator.
# See also CPackOptions.cmake.in, which gets modified at CMake time,
# and then used at CPack time.
############################################################################

set(CPACK_PACKAGE_NAME "${PYVINECOPULIB_PACKAGE_NAME}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CPACK_PACKAGE_NAME} - for doing something interesting.")
set(CPACK_PACKAGE_VENDOR "University College London.")
set(CPACK_PACKAGE_VERSION "${PYVINECOPULIB_VERSION_MAJOR}.${PYVINECOPULIB_VERSION_MINOR}.${PYVINECOPULIB_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_MAJOR "${PYVINECOPULIB_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PYVINECOPULIB_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PYVINECOPULIB_VERSION_PATCH}")
#set(CPACK_CREATE_DESKTOP_LINKS "QtVTKDemo")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_BINARY_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_BINARY_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_BINARY_DIR}/INSTALLATION.txt")
set(CPACK_PACKAGE_FILE_NAME "${PYVINECOPULIB_DEPLOY_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PYVINECOPULIB_DEPLOY_NAME}")
set(CPACK_MONOLITHIC_INSTALL ON)
