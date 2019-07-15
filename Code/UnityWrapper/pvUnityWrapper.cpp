/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "pvUnityWrapper.h"
#include "pvMyFunctions.h"

//-----------------------------------------------------------------------------
int AddTwoIntegers(int a, int b)
{
  return pv::MyFirstAddFunction(a, b);
}
