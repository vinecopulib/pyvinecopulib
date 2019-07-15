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

if(BUILD_QtVTKDemo AND VTK_FOUND AND Qt5_DIR)
  add_subdirectory(QtVTKDemo)
endif()

if(BUILD_QMLVTKDemo AND VTK_FOUND AND Qt5_DIR)
  add_subdirectory(QMLVTKDemo)
endif()

if(BUILD_QOpenGLDemo AND Qt5_DIR)
  add_subdirectory(QOpenGLDemo)
endif()
