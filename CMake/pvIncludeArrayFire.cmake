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

if(BUILD_ArrayFire)
  # Example of:
  #   (1) When SuperBuild builds ArrayFire it adds ArrayFire_DIR to CMAKE_PREFIX_PATH
  #       So, ArrayFire is found using ArrayFire's provided ArrayFireConfig.cmake.
  #       Its called 'config mode' when running find_package
  #       Its the preferred approach because ArrayFire can then control what is exposed.
  find_package(ArrayFire REQUIRED)
  include_directories(${ArrayFire_INCLUDE_DIRS})
  list(APPEND ALL_THIRD_PARTY_LIBRARIES ${ArrayFire_LIBRARIES})
  add_definitions(-DBUILD_ArrayFire)
  configure_file(${CMAKE_SOURCE_DIR}/Documentation/Licenses/ArrayFireAndTheirThirdPartyLicenses.txt ${CMAKE_BINARY_DIR}/LICENSE_ArrayFireAndTheirThirdPartyLicenses.txt)
  install(FILES ${CMAKE_BINARY_DIR}/LICENSE_ArrayFireAndTheirThirdPartyLicenses.txt DESTINATION . COMPONENT CONFIG)
endif()
