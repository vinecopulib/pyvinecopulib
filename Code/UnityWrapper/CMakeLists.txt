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

if(NOT BUILDING_GUIS)

  set(PYVINECOPULIB_UNITY_NAME ${PYVINECOPULIB_LIBRARY_NAME}Unity)

  set(PYVINECOPULIB_UNITY_SRCS
    pvUnityWrapper.cpp
  )

  set(PYVINECOPULIB_UNITY_HDRS
    pvUnityWrapper.h
  )

  # Thanks to: https://stackoverflow.com/questions/15120951/generating-a-bundle-file-with-cmake-on-mac-osx?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
  add_library(${PYVINECOPULIB_UNITY_NAME} MODULE ${PYVINECOPULIB_UNITY_HDRS} ${PYVINECOPULIB_UNITY_SRCS})
  if(APPLE)
    set_target_properties(${PYVINECOPULIB_UNITY_NAME} PROPERTIES BUNDLE TRUE)
  endif()

  target_link_libraries(${PYVINECOPULIB_UNITY_NAME} PRIVATE ${PYVINECOPULIB_LIBRARIES} ${ALL_THIRD_PARTY_LIBRARIES})
  install(TARGETS ${PYVINECOPULIB_UNITY_NAME} DESTINATION ${PYVINECOPULIB_INSTALL_LIB_DIR})

  if( NOT "${PYVINECOPULIB_UNITY_PLUGIN_DIR}" STREQUAL "")
    install(TARGETS ${PYVINECOPULIB_UNITY_NAME} DESTINATION ${PYVINECOPULIB_UNITY_PLUGIN_DIR})
  endif()

endif()
