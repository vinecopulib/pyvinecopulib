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

if(NOT BUILDING_GUIS)


  if (NOT PYVINECOPULIB_PYTHON_MODULE_NAME OR PYVINECOPULIB_PYTHON_MODULE_NAME STREQUAL "")
    set(PYVINECOPULIB_PYTHON_MODULE_NAME "" CACHE STRING "Basename for python module." FORCE)
  endif()

  set(PYVINECOPULIB_PYTHON_SRCS)

  if (BUILD_OpenCV)
    list(APPEND PYVINECOPULIB_PYTHON_SRCS pvLibPythonWithOpenCV.cpp) #Rename as appropriate.
    message("Including pyboostcvconverter, thanks to Gregory Kramida.")

    execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c
              "from __future__ import print_function\ntry: import numpy; print(numpy.get_include(), end='')\nexcept:pass\n"
              OUTPUT_VARIABLE NumPy_INCLUDE_DIR)
    execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c
              "from __future__ import print_function\ntry: import numpy; print(numpy.__version__, end='')\nexcept:pass\n"
      OUTPUT_VARIABLE NumPy_VERSION)
    message("Found NumPy_INCLUDE_DIR=${NumPy_INCLUDE_DIR}, NumPy_VERSION=${NumPy_VERSION}")
    include_directories(${NumPy_INCLUDE_DIR})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pyboostcvconverter/include)
    list(APPEND PYVINECOPULIB_PYTHON_SRCS pyboostcvconverter/src/pyboost_cv3_converter.cpp)
  else()
    list(APPEND PYVINECOPULIB_PYTHON_SRCS pvLibPython.cpp) # Rename as appropriate.
  endif()

  if (WIN32)
    if(NOT BUILD_SHARED)
      add_definitions(-DBOOST_PYTHON_STATIC_LIB)
      add_definitions(-DBOOST_NUMPY_STATIC_LIB)
    endif()
  endif()

  link_libraries(${ALL_LIBRARIES})

  PYTHON_ADD_MODULE(${PYVINECOPULIB_PYTHON_MODULE_NAME} ${PYVINECOPULIB_PYTHON_SRCS})
  if(APPLE)
    set_target_properties(${PYVINECOPULIB_PYTHON_MODULE_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif()

  install(TARGETS ${PYVINECOPULIB_PYTHON_MODULE_NAME} DESTINATION ${PYVINECOPULIB_INSTALL_LIB_DIR})

endif()
