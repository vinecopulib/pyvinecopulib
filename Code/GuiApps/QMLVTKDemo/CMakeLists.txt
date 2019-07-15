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

set(APP_NAME QMLVTKDemo)

set(APP_UIS
)

set(APP_SRCS
  Model/pvBaseModel.cpp
  Model/pvTriangleModel.cpp
  Model/pvCubeModel.cpp
  View/pvQOpenGLTriangleRenderer.cpp
  View/pvQQuickVTKView.cpp
)

set(APP_MOC_HDRS
  Model/pvBaseModel.h
  Model/pvTriangleModel.h
  Model/pvCubeModel.h
  View/pvQOpenGLTriangleRenderer.h
  View/pvQQuickVTKView.h
)

set(APP_HDRS
)

set(APP_RES
  Resources/QMLVTKDemo.rc
)

qt5_wrap_ui(QT5_UIS ${APP_UIS})
qt5_wrap_cpp(QT5_MOCS ${APP_MOC_HDRS})
qt5_add_resources(QT5_RES ${APP_RES})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Model)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/View)

# Install all Qt qml infrastructure.
# Note: Im trying to install everything under the bin/qml folder EXCEPT QtTest.
#       Im sure there's a nicer way of doing it using the QtQuickCompiler or something.
set(_qml_source_dir ${Qt5_DIR}/../../../qml)
set(_qml_dest_dir bin/qml)
if(APPLE)
  set(_qml_dest_dir ${APP_NAME}.app/Contents/MacOS/qml)
endif()
file(GLOB children RELATIVE ${_qml_source_dir} ${_qml_source_dir}/*)
foreach( c ${children})
  if( NOT ("QtTest" STREQUAL "${c}" ) )
    if( IS_DIRECTORY ${_qml_source_dir}/${c} )
      install(DIRECTORY ${_qml_source_dir}/${c} DESTINATION ${_qml_dest_dir})
    else()
      install(FILES ${_qml_source_dir}/${c} DESTINATION ${_qml_dest_dir})
    endif()
  endif()
endforeach()
qt5_add_resources(QT5_QML View/qml.qrc)


add_executable(${APP_NAME} MACOSX_BUNDLE
  ${APP_SRCS}
  ${APP_HDRS}
  ${QT5_UIS}
  ${QT5_MOCS}
  ${QT5_RES}
  ${QT5_QML}
  App/main.cpp
)

# Add libraries here if necessary.
set(APP_LIBS )
target_link_libraries(${APP_NAME} ${APP_LIBS} ${ALL_THIRD_PARTY_LIBRARIES} ${QT5_LINK_LIBRARIES})

# This sorts out the fixup_bundle.
pvMacroCreateGuiApplication(${APP_NAME} "${ADDITIONAL_SEARCH_PATHS}")
