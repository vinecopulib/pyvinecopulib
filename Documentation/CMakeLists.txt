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

file(GLOB DOXIN_FILES *.dox.in)

foreach(doxin ${DOXIN_FILES})
  get_filename_component(basename ${doxin} NAME_WE)
  configure_file(${doxin} ${CMAKE_BINARY_DIR}/Doxygen/${basename}.dox @ONLY)
endforeach()
