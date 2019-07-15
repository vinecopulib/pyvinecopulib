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

function(pvFunctionFixInstallNameForMac dir)
  file(GLOB dylibFiles ${dir}/lib/*.dylib)
  foreach(_dylib ${dylibFiles})
     message("Fixing install name for lib: ${_dylib}")
     execute_process(COMMAND install_name_tool -id ${_dylib} ${_dylib})
     foreach(_dep_dylib ${dylibFiles})
        get_filename_component(_dep_dylib_name ${_dep_dylib} NAME)
        execute_process(COMMAND install_name_tool -change ${_dep_dylib_name} \@loader_path/${_dep_dylib_name} ${_dylib})
     endforeach()
  endforeach()
endfunction()



