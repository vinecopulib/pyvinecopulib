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

find_path(Eigen_INCLUDE_DIR
  NAMES Eigen/Eigen
  PATHS ${Eigen_DIR} ${CMAKE_PREFIX_PATH}
  PATH_SUFFIXES include include/eigen3
)

if (NOT TARGET Eigen)
  add_library(Eigen INTERFACE IMPORTED GLOBAL)
  set_property(TARGET Eigen APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Eigen_INCLUDE_DIR})
endif()

find_package_handle_standard_args(Eigen
  FOUND_VAR Eigen_FOUND
  REQUIRED_VARS Eigen_INCLUDE_DIR
)

