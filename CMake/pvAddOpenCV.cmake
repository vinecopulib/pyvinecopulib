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

option(BUILD_OpenCV "Build OpenCV." OFF)

if(BUILD_VTK AND APPLE AND BUILD_OpenCV AND "${VTK_VERSION}" STREQUAL "${DEFAULT_VTK_VERSION}")
  set(VTK_VERSION "${FALLBACK_VTK_VERSION}")
  message("Forcing VTK_VERSION to ${VTK_VERSION} as you are on Mac OSX and BUILD_OpenCV is on.")
endif()
