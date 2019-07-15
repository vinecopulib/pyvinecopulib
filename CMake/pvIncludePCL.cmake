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

if(BUILD_PCL)
  # Example of:
  #   (1) When SuperBuild builds PCL it adds PCL_DIR to CMAKE_PREFIX_PATH
  #       So, PCL is found using PCL's provided PCLConfig.cmake.
  #       Its called 'config mode' when running find_package
  #       Its the preferred approach because PCL can then control what is exposed.
  set(PCL_FIND_QUIETLY ON)
  find_package(PCL REQUIRED)
  include_directories(${PCL_INCLUDE_DIRS})
  link_directories(${PCL_LIBRARY_DIRS})
  add_definitions(${PCL_DEFINITIONS})
  list(APPEND ALL_THIRD_PARTY_LIBRARIES ${PCL_LIBRARIES})

  # This appears to be missing from the list of PCL_LIBRARIES,
  # but I don't yet know why.
  find_library(_pcl_io_ply_LIBRARY
    NAMES pcl_io_ply
    PATHS ${PCL_LIBRARY_DIRS}
    NO_DEFAULT_PATH
  )
  if(_pcl_io_ply_LIBRARY)
    list(APPEND ALL_THIRD_PARTY_LIBRARIES ${_pcl_io_ply_LIBRARY})
  endif()

  add_definitions(-DBUILD_PCL)
  list(APPEND ADDITIONAL_SEARCH_PATHS "${PCL_LIBRARY_DIRS}/../${_library_sub_dir}")
  list(APPEND ADDITIONAL_SEARCH_PATHS "${FLANN_DIR}/${_library_sub_dir}")
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Licenses/PCL.txt ${CMAKE_BINARY_DIR}/LICENSE_PCL.txt)
  install(FILES ${CMAKE_BINARY_DIR}/LICENSE_PCL.txt DESTINATION . COMPONENT CONFIG)
endif()
