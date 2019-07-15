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

configure_file(${CMAKE_SOURCE_DIR}/UsePyVinecopulib.cmake.in ${CMAKE_BINARY_DIR}/UsePyVinecopulib.cmake @ONLY IMMEDIATE)
configure_file(${CMAKE_SOURCE_DIR}/PyVinecopulibConfig.cmake.in ${CMAKE_BINARY_DIR}/PyVinecopulibConfig.cmake @ONLY IMMEDIATE)
if(NOT BUILDING_GUIS)
  install(FILES ${CMAKE_BINARY_DIR}/UsePyVinecopulib.cmake DESTINATION . COMPONENT CONFIG)
  install(FILES ${CMAKE_BINARY_DIR}/PyVinecopulibConfig.cmake DESTINATION . COMPONENT CONFIG)
endif()
