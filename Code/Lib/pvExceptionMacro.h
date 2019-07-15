/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvExceptionMacro_h
#define pvExceptionMacro_h

#include "pvException.h"

#define pvExceptionThrow() throw pv::Exception(__FILE__,__LINE__)

#endif
