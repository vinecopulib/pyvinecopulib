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

option(BUILD_Docs "Build Docs using Doxygen." OFF)

# See: https://www.stack.nl/~dimitri/doxygen/manual/external.html
# When you do find_package on an external project, you should also append
# ext1/ext1.tag=../../ext1/html i.e. name=path
# to this variable:
set(PYVINECOPULIB_EXTERNAL_DOXYGEN_TAGFILES)
