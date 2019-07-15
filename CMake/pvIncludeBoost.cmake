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

if(BUILD_Boost OR BUILD_Python_Boost)
  # Example of:
  #   (1) Standard, widely used library
  #   (2) Mostly header library, some compiled libraries
  #   (3) Your CMake probably has a standard FindBoost.cmake in its distribution,
  #       so Boost is probably found by CMake's own FindBoost.cmake.
  #       which is an example of the 'module' mode.
  set(BOOST_LIBRARYDIR ${BOOST_ROOT}/lib)
  set(BOOST_LIBRARYDIR ${BOOST_ROOT}/lib)
  set(Boost_LIBRARY_DIR_DEBUG ${BOOST_ROOT}/lib)
  set(Boost_LIBRARY_DIR_RELEASE ${BOOST_ROOT}/lib)
  set(BOOST_INCLUDEDIR ${BOOST_ROOT}/include)
  #set(Boost_DEBUG ON)
  set(Boost_NO_SYSTEM_PATHS ON)  # Notice: This enables us to turn off System Paths
  set(Boost_NO_BOOST_CMAKE ON)   # Notice: We can tell CMake not to assume this is the CMake-ified version of the boost project.
  if(BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS OFF)
    set(Boost_USE_STATIC_RUNTIME OFF)
  else()
    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_STATIC_RUNTIME ON)
  endif()
  find_package(Boost 1.67 EXACT COMPONENTS ${PYVINECOPULIB_BOOST_LIBS} REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
  message("Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
  message("Boost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}")
  if(WIN32)
    if(WITHIN_SUBBUILD)
      add_definitions(-DBoost_LIB_DIAGNOSTIC_DEFINITIONS)  # To get debug messages
      add_definitions(-DBOOST_ALL_NO_LIB)                  # To stop auto-linking, which seems to be adding "lib" as library prefix in .obj files.
    endif()
    if(BUILD_SHARED)
      list(APPEND ALL_COMPILE_OPTIONS -DBOOST_ALL_DYN_LINK)
    endif()
  endif()
  list(APPEND ALL_THIRD_PARTY_LIBRARIES ${Boost_LIBRARIES})
  add_definitions(-DBUILD_Boost)
  add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
  list(APPEND ADDITIONAL_SEARCH_PATHS "${BOOST_ROOT}/${_library_sub_dir}")
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Licenses/Boost.txt ${CMAKE_BINARY_DIR}/LICENSE_Boost.txt)
  install(FILES ${CMAKE_BINARY_DIR}/LICENSE_Boost.txt DESTINATION . COMPONENT CONFIG)

  if(BUILD_SHARED AND WITHIN_SUBBUILD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_LIB_PREFIX=\"\"")
  endif()
  
endif()
