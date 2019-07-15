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

if(PYVINECOPULIB_USE_OPENMP)
  # Borrowed from PCL-1.8, and modified
  find_package(OpenMP)

  # Note: Convention would say that for `find_package(OpenMP)`, the variable should be OpenMP_FOUND.
  # Also: The top of the /usr/local/share/cmake/Modules/FindOpenMP.cmake says it should be OpenMP_FOUND.
  # But : Both PCL and OpenCV use OPENMP_FOUND, then I checked FindOpenMP.cmake and its OPENMP_FOUND.

  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    list(APPEND ALL_THIRD_PARTY_LIBRARIES ${OpenMP_CXX_LIBRARIES})
    message("OpenMP found, CMAKE_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
    message("OpenMP found, OpenMP_CXX_LIBRARIES=${OpenMP_CXX_LIBRARIES}")
    if(MSVC)
      if(MSVC_VERSION EQUAL 1500)
        set(OPENMP_DLL VCOMP90)
      elseif(MSVC_VERSION EQUAL 1600)
        set(OPENMP_DLL VCOMP100)
      elseif(MSVC_VERSION EQUAL 1700)
        set(OPENMP_DLL VCOMP110)
      elseif(MSVC_VERSION EQUAL 1800)
        set(OPENMP_DLL VCOMP120)
      elseif(MSVC_VERSION EQUAL 1900)
        set(OPENMP_DLL VCOMP140)
      elseif(MSVC_VERSION EQUAL 1910)
        set(OPENMP_DLL VCOMP140)
      endif()
      if(OPENMP_DLL)
        set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DELAYLOAD:${OPENMP_DLL}D.dll")
        set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DELAYLOAD:${OPENMP_DLL}.dll")
      else(OPENMP_DLL)
        message(WARNING "Delay loading flag for OpenMP DLL is invalid.")
      endif(OPENMP_DLL)
    endif(MSVC)
  endif()
endif()
