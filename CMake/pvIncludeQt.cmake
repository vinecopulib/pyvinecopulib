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

if(PYVINECOPULIB_USE_QT)
  set(_qt_components Core Concurrent PrintSupport Script Sql Svg Xml XmlPatterns)
  if(BUILD_QtVTKDemo OR BUILD_QOpenGLDemo)
    list(APPEND _qt_components OpenGL Gui Widgets UiTools Help)
  endif()
  if(BUILD_QMLVTKDemo)
    list(APPEND _qt_components OpenGL Gui Quick Qml)
  endif()
  if(UNIX AND NOT APPLE)
    if(BUILD_SHARED_LIBS)
      list(APPEND _qt_components X11Extras)
    endif()
  endif()
  find_package(Qt5 QUIET COMPONENTS ${_qt_components})
  if(Qt5_DIR)
    message(STATUS "Found Qt5: ${Qt5_DIR}")
    if (Qt5Core_VERSION VERSION_LESS 5.5.0 AND "${VTK_VERSION}" STREQUAL "${DEFAULT_VTK_VERSION}")
      message("Forcing VTK_BACKEND to OpenGL due to Qt version < 5.5.0")
      set(VTK_BACKEND "OpenGL")
    endif()
    get_filename_component(_Qt5_DIR "${Qt5_DIR}/../../../" ABSOLUTE)
    list(FIND CMAKE_PREFIX_PATH "${_Qt5_DIR}" _result)
    if(_result LESS 0)
      set(CMAKE_PREFIX_PATH "${_Qt5_DIR};${CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE)
    endif()
    set(PYVINECOPULIB_PREFIX_PATH ${_Qt5_DIR})
    foreach(_component ${_qt_components})
      find_package(Qt5${_component} REQUIRED QUIET)
      include_directories(${Qt5${_component}_INCLUDE_DIRS})
      add_definitions(${Qt5${_component}_DEFINITIONS})
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5${_component}_EXECUTABLE_COMPILE_FLAGS}")
      list(APPEND QT5_LINK_LIBRARIES Qt5::${_component})
    endforeach()
    set(_qmake_location )
    if(TARGET ${Qt5Core_QMAKE_EXECUTABLE})
      get_property(_qmake_location TARGET ${Qt5Core_QMAKE_EXECUTABLE}
                   PROPERTY IMPORT_LOCATION)
    endif()
    if(_qmake_location)
      if(NOT _qt_install_libs)
        if(WIN32)
          execute_process(COMMAND ${_qmake_location} -query QT_INSTALL_BINS
                          OUTPUT_VARIABLE _qt_install_libs
                          OUTPUT_STRIP_TRAILING_WHITESPACE)
        else()
          execute_process(COMMAND ${_qmake_location} -query QT_INSTALL_LIBS
                          OUTPUT_VARIABLE _qt_install_libs
                          OUTPUT_STRIP_TRAILING_WHITESPACE)
        endif()
        file(TO_CMAKE_PATH "${_qt_install_libs}" _qt_install_libs)
        set(_qt_install_libs ${_qt_install_libs} CACHE INTERNAL "Qt library installation prefix" FORCE)
      endif()
      if(_qt_install_libs)
        list(APPEND ADDITIONAL_SEARCH_PATHS ${_qt_install_libs})
      endif()
    else()
      message(WARNING "The qmake executable could not be found.")
    endif()
  endif()
endif()

if(PYVINECOPULIB_USE_QT AND NOT Qt5_DIR)
  message(FATAL_ERROR "Qt was required but not found")
endif()
