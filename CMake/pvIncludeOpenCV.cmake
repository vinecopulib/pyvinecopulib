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

if(BUILD_OpenCV)
  # Example of:
  #   (1) When SuperBuild builds OpenCV it adds OpenCV_DIR to CMAKE_PREFIX_PATH
  #       So, OpenCV is found using OpenCV's provided OpenCVConfig.cmake.
  #       Its called 'config mode' when running find_package
  #       Its the preferred approach because OpenCV can then control what is exposed.
  find_package(OpenCV REQUIRED)
  include_directories(${OpenCV_INCLUDE_DIRS})
  list(APPEND ALL_THIRD_PARTY_LIBRARIES ${OpenCV_LIBS})
  add_definitions(-DBUILD_OpenCV)
  if(WIN32)
    list(APPEND ADDITIONAL_SEARCH_PATHS "${OpenCV_LIB_PATH}/../${_library_sub_dir}")
  else()
    list(APPEND ADDITIONAL_SEARCH_PATHS "${OpenCV_INSTALL_PATH}/${_library_sub_dir}")
  endif()
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Licenses/OpenCV.txt ${CMAKE_BINARY_DIR}/LICENSE_OpenCV.txt)
  install(FILES ${CMAKE_BINARY_DIR}/LICENSE_OpenCV.txt DESTINATION . COMPONENT CONFIG)
endif()
