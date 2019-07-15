/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/
/**
* \file pvUnityWrapper.h
* \brief Wrapper for Unity3D to load simple functionality as plugin.
* \ingroup utilities
*
* Note: this example is based on:
* https://docs.unity3d.com/Manual/PluginsForDesktop.html
*
* For a more advanced plugin see:
* https://docs.unity3d.com/Manual/NativePluginInterface.html
*/

/**
 * Unity plugin's must have C-style linkage to avoid name mangling.
 */
extern "C" {

/**
 * \brief C-style wrapper for pv::MyFirstAddFunction
 */
int AddTwoIntegers(int a, int b);

}
