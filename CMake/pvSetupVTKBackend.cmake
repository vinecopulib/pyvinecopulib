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

if("${VTK_VERSION}" STREQUAL "${FALLBACK_VTK_VERSION}" AND "${VTK_BACKEND}" STREQUAL "${DEFAULT_VTK_BACKEND}")
  set(VTK_BACKEND "OpenGL")
  message("Forcing VTK_BACKEND to OpenGL instead of OpenGL2, due to VTK_VERSION=${VTK_VERSION}")
endif()
