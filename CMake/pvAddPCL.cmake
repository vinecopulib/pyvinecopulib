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

option(BUILD_PCL "Build PCL." OFF)
if(BUILD_PCL)
  find_package(OpenGL REQUIRED)
  option(BUILD_PCL_VIS "Build PCL Visualisation tools." OFF)
endif()

if(BUILD_PCL AND NOT BUILD_Boost)
  set(BUILD_Boost ON CACHE BOOL "Build Boost." FORCE)
  message("Forcing BUILD_Boost to ON due to BUILD_PCL being ON.")
endif()

if(BUILD_PCL AND NOT BUILD_Eigen)
  set(BUILD_Eigen ON CACHE BOOL "Build Eigen." FORCE)
  message("Forcing BUILD_Eigen to ON due to BUILD_PCL being ON.")
endif()

if(BUILD_PCL AND NOT BUILD_FLANN)
  set(BUILD_FLANN ON CACHE BOOL "Build FLANN." FORCE)
  message("Forcing BUILD_FLANN to ON due to BUILD_PCL being ON.")
endif()

if(BUILD_PCL_VIS AND NOT BUILD_PCL)
  set(BUILD_PCL_VIS OFF CACHE BOOL "Build PCL Visualisation tools." FORCE)
  message("Forcing BUILD_PCL_VIS to OFF due to BUILD_PCL being OFF.")
endif()

if(BUILD_PCL_VIS AND BUILD_Python_Boost)
  set(BUILD_PCL_VIS OFF CACHE BOOL "Build PCL Visualisation tools." FORCE)
  message("Forcing BUILD_PCL_VIS to OFF due to BUILD_Python_Boost being ON.")
endif()

if(BUILD_PCL_VIS AND BUILD_Python_PyBind)
  set(BUILD_PCL_VIS OFF CACHE BOOL "Build PCL Visualisation tools." FORCE)
  message("Forcing BUILD_PCL_VIS to OFF due to BUILD_Python_PyBind being ON.")
endif()

# due to https://github.com/PointCloudLibrary/pcl/issues/712
if(BUILD_VTK AND APPLE AND BUILD_PCL AND BUILD_PCL_VIS AND "${VTK_VERSION}" STREQUAL "${DEFAULT_VTK_VERSION}")
  set(VTK_VERSION "${FALLBACK_VTK_VERSION}")
  message("Forcing VTK_VERSION to ${VTK_VERSION} as you are on Mac OSX and both BUILD_PCL and BUILD_PCL_VIS are on.")
endif()
