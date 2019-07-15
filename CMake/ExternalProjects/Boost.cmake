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

set(PYVINECOPULIB_USE_Boost 1)
set(PYVINECOPULIB_USE_Boost_LIBRARIES ${PYVINECOPULIB_BOOST_LIBS})
set(PYVINECOPULIB_BUILD_PYTHON_BINDINGS OFF)

#-----------------------------------------------------------------------------
# Boost
#-----------------------------------------------------------------------------
if(NOT BUILD_Boost AND NOT BUILD_Python_Boost)
  return()
endif()

# Sanity checks
if(DEFINED BOOST_ROOT AND NOT EXISTS ${BOOST_ROOT})
  message(FATAL_ERROR "BOOST_ROOT variable is defined but corresponds to non-existing directory")
endif()

set(version "1_67_0")
set(location "${NIFTK_EP_TARBALL_LOCATION}/boost_${version}.tar.gz")
mpMacroDefineExternalProjectVariables(Boost ${version} ${location})
set(proj_DEPENDENCIES )

string(REPLACE "^^" ";" PYVINECOPULIB_USE_Boost_LIBRARIES "${PYVINECOPULIB_USE_Boost_LIBRARIES}")

if(NOT DEFINED BOOST_ROOT AND NOT PYVINECOPULIB_USE_SYSTEM_Boost)

  set(_boost_version 1_67)
  set(_boost_install_include_dir include/boost)
  if(WIN32)
    set(_boost_install_include_dir include/boost-${_boost_version}/boost)
  endif()

  set(_boost_libs )
  set(_with_boost_libs )
  set(_install_lib_dir )

  # Set the boost root to the libraries install directory
  set(BOOST_ROOT "${proj_INSTALL}")

  if(PYVINECOPULIB_USE_Boost_LIBRARIES)
    string(REPLACE ";" "," _boost_libs "${PYVINECOPULIB_USE_Boost_LIBRARIES}")
    foreach(_boost_lib ${PYVINECOPULIB_USE_Boost_LIBRARIES})
      list(APPEND _with_boost_libs --with-${_boost_lib})
	  if("${_boost_lib}" STREQUAL "python" OR "${_boost_lib}" STREQUAL "python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
	    set(PYVINECOPULIB_BUILD_PYTHON_BINDINGS ON)
	  endif()
    endforeach()
  endif()

  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(_boost_address_model "address-model=64")
  else()
    set(_boost_address_model "address-model=32")
  endif()

  if(WIN32)
    set(_shell_extension .bat)
    set(_boost_layout)
    if(MSVC)
      if(MSVC_VERSION EQUAL 1600)
        set(_boost_with_toolset "vc10")
        set(_boost_toolset "msvc-10.0")
      elseif(MSVC_VERSION EQUAL 1700)
        set(_boost_with_toolset "vc11")
        set(_boost_toolset "msvc-11.0")
      elseif(MSVC_VERSION EQUAL 1800)
        set(_boost_with_toolset "vc12")
        set(_boost_toolset "msvc-12.0")
      elseif(MSVC_VERSION EQUAL 1900)
        set(_boost_with_toolset "vc14")
        set(_boost_toolset "msvc-14.0")
      endif()
    endif()
    set(_install_lib_dir "--libdir=<INSTALL_DIR>/bin")
    set(WIN32_CMAKE_SCRIPT ${proj_CONFIG}/MoveBoostLibsToLibDirForWindows.cmake)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMake/ExternalProjects/MoveBoostLibsToLibDirForWindows.cmake.in ${WIN32_CMAKE_SCRIPT} @ONLY)
    set(_windows_move_libs_cmd COMMAND ${CMAKE_COMMAND} -P ${WIN32_CMAKE_SCRIPT})
  else()
    set(_shell_extension .sh)
    set(_boost_layout "--layout=tagged")
  endif()

  if(UNIX AND NOT APPLE)
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(_boost_with_toolset "gcc")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(_boost_with_toolset "clang")
    else()
      message(FATAL_ERROR "Compiler '${CMAKE_CXX_COMPILER_ID}' not supported. Use GNU or Clang instead.")
    endif()
    get_filename_component(_cxx_compiler_name "${CMAKE_CXX_COMPILER}" NAME)
    string(REGEX MATCH "^[0-9]+\\.[0-9]+" _compiler_version "${CMAKE_CXX_COMPILER_VERSION}")
    if(_cxx_compiler_name MATCHES "${_compiler_version}")
      set(_boost_toolset "${_boost_with_toolset}-${_compiler_version}")
    endif()
  endif()

  if(_boost_toolset)
    set(_boost_toolset "--toolset=${_boost_toolset}")
  endif()

  set(APPLE_CLANG_FLAGS)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" AND "${CMAKE_OSX_DEPLOYMENT_TARGET}" STREQUAL "10.8" )
    set(APPLE_CLANG_FLAGS toolset=clang cxxflags="-stdlib=libstdc++" linkflags="-stdlib=libstdc++")
  endif()

  set (APPLE_SYSROOT_FLAG)
  if(APPLE)
    set(APPLE_CMAKE_SCRIPT ${proj_CONFIG}/ChangeBoostLibsInstallNameForMac.cmake)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMake/ExternalProjects/ChangeBoostLibsInstallNameForMac.cmake.in ${APPLE_CMAKE_SCRIPT} @ONLY)
    set(_macos_change_install_name_cmd COMMAND ${CMAKE_COMMAND} -P ${APPLE_CMAKE_SCRIPT})

    # Set OSX_SYSROOT
    if (NOT ${CMAKE_OSX_SYSROOT} STREQUAL "")
      set (APPLE_SYSROOT_FLAG --sysroot=${CMAKE_OSX_SYSROOT})
    endif()
  endif()

  if(PYVINECOPULIB_BUILD_PYTHON_BINDINGS)
    set(_python_version "--with-python-version=${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
    get_filename_component(_python_root_dir "${PYTHON_EXECUTABLE}" DIRECTORY)
    set(_python_root "--with-python-root=${_python_root_dir}")
    message("Building boost::python with:${_python_version};${_python_root}")
  endif()

  set(_boost_variant "$<$<CONFIG:Debug>:debug>$<$<CONFIG:Release>:release>")
  set(_boost_link shared)
  if(NOT BUILD_SHARED_LIBS)
    set(_boost_link static)
  endif()
  set(_boost_cxxflags )
  if(CMAKE_CXX_FLAGS OR PYVINECOPULIB_CXX11_FLAG OR PythonLibs_FOUND)
    set(_boost_cxxflags "cxxflags=${PYVINECOPULIB_CXX11_FLAG} ${CMAKE_CXX_FLAGS}")
    if(WIN32)
      set(_boost_cxxflags "${_boost_cxxflags} /I${PYTHON_INCLUDE_DIRS}")
    else()
      set(_boost_cxxflags "${_boost_cxxflags} -I${PYTHON_INCLUDE_DIRS}")
    endif()
  endif()
  set(_boost_linkflags )
  if(BUILD_SHARED_LIBS AND _install_rpath_linkflag)
    set(_boost_linkflags "linkflags=${_install_rpath_linkflag}")
  endif()

  set(_build_cmd "<SOURCE_DIR>/b2"
      ${APPLE_CLANG_FLAGS}
      ${APPLE_SYSROOT_FLAG}
      ${_boost_toolset}
      ${_boost_layout}
      "--prefix=<INSTALL_DIR>"
      ${_install_lib_dir}
      ${_with_boost_libs}
      # Use the option below to view the shell commands (for debugging)
      #-d+4
      -d 0
      variant=${_boost_variant}
      link=${_boost_link}
      ${_boost_cxxflags}
      ${_boost_linkflags}
      ${_boost_address_model}
      threading=multi
      runtime-link=${_boost_link}
      --ignore-site-config
      -q
  )

  if(PYVINECOPULIB_USE_Boost_LIBRARIES)
    set(_boost_build_cmd BUILD_COMMAND ${_build_cmd})
    set(_install_cmd ${_build_cmd} install ${_macos_change_install_name_cmd} ${_windows_move_libs_cmd})
  else()
    set(_boost_build_cmd BUILD_COMMAND ${CMAKE_COMMAND} -E echo "no binary libraries")
    set(_install_cmd ${CMAKE_COMMAND} -E echo "copying Boost header..."
           COMMAND ${CMAKE_COMMAND} -E copy_directory "<SOURCE_DIR>/boost" "<INSTALL_DIR>/${_boost_install_include_dir}")
  endif()

  ExternalProject_Add(${proj}
    LIST_SEPARATOR ^^
    PREFIX ${proj_CONFIG}
    SOURCE_DIR ${proj_SOURCE}
    # Boost needs in-source builds
    BINARY_DIR ${proj_SOURCE}
    INSTALL_DIR ${proj_INSTALL}
    URL ${proj_LOCATION}
    URL_MD5 ${proj_CHECKSUM}
    CONFIGURE_COMMAND "<SOURCE_DIR>/bootstrap${_shell_extension}"
      --with-toolset=${_boost_with_toolset}
      --with-libraries=${_boost_libs}
      ${_python_version}
      "${_python_root}"
      "--prefix=<INSTALL_DIR>"
    ${_boost_build_cmd}
    INSTALL_COMMAND ${_install_cmd}
    DEPENDS ${proj_DEPENDENCIES}
    )

  set(PYVINECOPULIB_PREFIX_PATH ${proj_INSTALL}^^${PYVINECOPULIB_PREFIX_PATH})

  ExternalProject_Get_Property(${proj} install_dir)

  if(WIN32)
    set(BOOST_LIBRARYDIR "${install_dir}/lib")
  endif()

  # Manual install commands (for a PYVINECOPULIB super-build install)
  # until the Boost CMake system is used.

  # We just copy the include directory
  install(DIRECTORY "${install_dir}/${_boost_install_include_dir}"
          DESTINATION "include"
          COMPONENT dev
         )

  if(PYVINECOPULIB_USE_Boost_LIBRARIES)
    # Copy the boost libraries
    file(GLOB _boost_libs
         "${install_dir}/lib/libboost*.so*"
         "${install_dir}/lib/libboost*.dylib")
    install(FILES ${_boost_libs}
            DESTINATION "lib"
            COMPONENT runtime)

    file(GLOB _boost_libs
         "${install_dir}/bin/libboost*.dll")
    install(FILES ${_boost_libs}
            DESTINATION "bin"
            COMPONENT runtime)

    file(GLOB _boost_libs
         "${install_dir}/lib/libboost*.lib"
         "${install_dir}/lib/libboost*.a")
    install(FILES ${_boost_libs}
            DESTINATION "lib"
            COMPONENT dev)
  endif()

  message("SuperBuild loading Boost from ${proj_INSTALL}")

else()

  mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

endif()
