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

add_subdirectory(Model)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Model)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/Model)

add_subdirectory(View)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/View)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/View)

add_subdirectory(Controller)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Controller)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/Controller)

add_subdirectory(App)
