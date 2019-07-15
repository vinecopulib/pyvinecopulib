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

if(BUILD_glog)

  # This is generated into build folder, so we find it before any system ones.
  configure_file(${CMAKE_SOURCE_DIR}/CMake/Findglog.cmake.in ${CMAKE_BINARY_DIR}/Findglog.cmake @ONLY)

  # Example of: Small library where Findglog.cmake is
  #             generated into the build folder to always pick up our one.
  find_package(glog REQUIRED)
  include_directories(${glog_INCLUDE_DIR})
  list(APPEND ALL_THIRD_PARTY_LIBRARIES ${glog_LIBRARY})
  add_definitions(-DBUILD_glog)
  list(APPEND ADDITIONAL_SEARCH_PATHS "${glog_DIR}/${_library_sub_dir}")
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Licenses/glog.txt ${CMAKE_BINARY_DIR}/LICENSE_glog.txt)
  install(FILES ${CMAKE_BINARY_DIR}/LICENSE_glog.txt DESTINATION . COMPONENT CONFIG)
  
endif()
