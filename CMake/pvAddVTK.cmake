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

option(BUILD_VTK "Build VTK." OFF)

set(DEFAULT_VTK_VERSION "v8.2.0")
set(FALLBACK_VTK_VERSION "v6.1.0") # Should only be needed on Mac, if PCL visualization tools are on, or OpenCV visualisation on.
set(DEFAULT_VTK_BACKEND "OpenGL2")
set(VTK_VERSION "${DEFAULT_VTK_VERSION}")
set(VTK_BACKEND "${DEFAULT_VTK_BACKEND}")
