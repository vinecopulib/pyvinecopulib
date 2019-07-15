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

if(WIN32)
  set(VS_SOLUTION_FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.sln")
  foreach(VS_BUILD_TYPE ${CMAKE_CONFIGURATION_TYPES})
    configure_file("${CMAKE_SOURCE_DIR}/CMake/StartVS.bat.in" ${PROJECT_BINARY_DIR}/StartVS_${VS_BUILD_TYPE}.bat @ONLY)
    message( "CreateWindowsBatchScript: Creating ${PROJECT_BINARY_DIR}/StartVS_${VS_BUILD_TYPE}.bat" )
  endforeach()
endif(WIN32)
