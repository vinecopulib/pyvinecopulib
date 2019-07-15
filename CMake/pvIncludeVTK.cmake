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

if(BUILD_VTK)
  # Example of:
  #   (1) Large Rendering library - oooooh, pretty pictures etc.
  #   (2) When SuperBuild builds VTK it adds VTK_DIR to CMAKE_PREFIX_PATH
  #       So, VTK is found using VTK's provided VTKConfig.cmake.
  #       Its called 'config mode' when running find_package
  #       Its the preferred approach because VTK can then control what is exposed.
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
  list(APPEND ALL_THIRD_PARTY_LIBRARIES ${VTK_LIBRARIES})
  add_definitions(-DBUILD_VTK)
  add_definitions(-DBUILD_VTK_${VTK_BACKEND})
  list(APPEND ADDITIONAL_SEARCH_PATHS "${VTK_INSTALL_PREFIX}/${_library_sub_dir}")
  if(UNIX AND NOT APPLE)
    list(APPEND ADDITIONAL_SEARCH_PATHS "${VTK_INSTALL_PREFIX}/lib64")
  endif()
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Licenses/VTK.txt ${CMAKE_BINARY_DIR}/LICENSE_VTK.txt)
  install(FILES ${CMAKE_BINARY_DIR}/LICENSE_VTK.txt DESTINATION . COMPONENT CONFIG)
endif()
