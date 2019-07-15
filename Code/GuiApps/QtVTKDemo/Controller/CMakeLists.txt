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

#########################
# Put .ui files here.
#########################
set(QTVTKCONTROLLER_UIS
  pvTwinSliderWidget.ui
  pvControlPanelWidget.ui
)


#########################
# Put .cpp files here.
#########################
set(QTVTKCONTROLLER_SRCS
  pvTwinSliderWidget.cpp
  pvControlPanelWidget.cpp
)


###################################
# Put .h files that need moc here.
###################################
SET(QTVTKCONTROLLER_MOC_HDRS
  pvTwinSliderWidget.h
  pvControlPanelWidget.h
)


##########################################
# Put .h files that do not need moc here.
##########################################
SET(QTVTKCONTROLLER_HDRS
  pvQtVTKControllerWin32ExportHeader.h
)


#########################
# Resource file.
#########################
SET(QTVTKCONTROLLER_RES Resources/qtvtkcontroller.qrc)


################################
# Qt commands provided by CMake.
################################
qt5_wrap_ui(QT5_UIS ${QTVTKCONTROLLER_UIS} )
qt5_wrap_cpp(QT5_MOCS ${QTVTKCONTROLLER_MOC_HDRS} )
qt5_add_resources(QT5_RES ${QTVTKCONTROLLER_RES})

################################
# Bundle it all together.
################################
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(qtvtkcontroller
  ${QTVTKCONTROLLER_SRCS}
  ${QTVTKCONTROLLER_HDRS}
  ${QT5_UIS}
  ${QT5_MOCS}
  ${QT5_RES}
)


#################################
# Stuff for windows shared libs.
#################################
if(BUILD_SHARED_LIBS)
  if(WIN32)
    add_definitions(-DPYVINECOPULIB_QTVTKCONTROLLER_WINDOWS_EXPORT)
  endif(WIN32)
  PYVINECOPULIB_INSTALL_LIBRARY(qtvtkcontroller)
endif(BUILD_SHARED_LIBS)


#################################
# Link the library.
#################################
target_link_libraries(qtvtkcontroller PRIVATE ${QT5_LINK_LIBRARIES} ${ALL_THIRD_PARTY_LIBRARIES})

