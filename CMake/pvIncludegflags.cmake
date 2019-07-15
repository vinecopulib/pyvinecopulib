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

if(BUILD_gflags)

  # This is generated into build folder, so we find it before any system ones.
  configure_file(${CMAKE_SOURCE_DIR}/CMake/Findgflags.cmake.in ${CMAKE_BINARY_DIR}/Findgflags.cmake @ONLY)

  # Example of: Small library where Findgflags.cmake is
  #             generated into the build folder to always pick up our one.
  find_package(gflags REQUIRED)
  find_package(Threads REQUIRED)
  include_directories(${gflags_INCLUDE_DIR})
  list(APPEND ALL_THIRD_PARTY_LIBRARIES ${gflags_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
  add_definitions(-DBUILD_gflags)
  if(WIN32)
    add_definitions(-DGOOGLE_GLOG_DLL_DECL=)
    list(APPEND ALL_THIRD_PARTY_LIBRARIES Shlwapi)
  endif()
  if(WIN32)
    list(APPEND ADDITIONAL_SEARCH_PATHS "${gflags_DIR}/Lib")
  else()
    list(APPEND ADDITIONAL_SEARCH_PATHS "${gflags_DIR}/${_library_sub_dir}")
  endif()
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Licenses/gflags.txt ${CMAKE_BINARY_DIR}/LICENSE_gflags.txt)
  install(FILES ${CMAKE_BINARY_DIR}/LICENSE_gflags.txt DESTINATION . COMPONENT CONFIG)
endif()
