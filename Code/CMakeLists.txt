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

add_subdirectory(Lib)
add_subdirectory(CommandLineApps) # Must go before GuiApps for fixup_bundle to work.
add_subdirectory(GuiApps)

if(BUILD_Python_Boost)
  add_subdirectory(PythonBoost)
endif()

if(BUILD_Python_PyBind)
  add_subdirectory(PythonPyBind)
endif()

if(BUILD_UNITY_WRAPPER)
  add_subdirectory(UnityWrapper)
endif()