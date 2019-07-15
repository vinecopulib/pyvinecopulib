/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvWin32ExportHeader_h
#define pvWin32ExportHeader_h

/**
* \file pvWin32ExportHeader.h
* \brief Header to sort Windows dllexport/dllimport.
*/

#if (defined(_WIN32) || defined(WIN32)) && !defined(PYVINECOPULIB_STATIC)
  #ifdef PYVINECOPULIB_WINDOWS_EXPORT
    #define PYVINECOPULIB_WINEXPORT __declspec(dllexport)
  #else
    #define PYVINECOPULIB_WINEXPORT __declspec(dllimport)
  #endif  /* PYVINECOPULIB_WINEXPORT */
#else
/* linux/mac needs nothing */
  #define PYVINECOPULIB_WINEXPORT
#endif

#endif
