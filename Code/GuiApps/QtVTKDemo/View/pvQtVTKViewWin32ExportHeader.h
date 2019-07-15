/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvQtVTKViewWin32ExportHeader_h
#define pvQtVTKViewWin32ExportHeader_h

/**
* \file pvQtVTKViewWin32ExportHeader.h
* \brief Header to sort Windows dllexport/dllimport.
*/

#if (defined(_WIN32) || defined(WIN32)) && !defined(PYVINECOPULIB_STATIC)
  #ifdef PYVINECOPULIB_QTVTKVIEW_WINDOWS_EXPORT
    #define PYVINECOPULIB_QTVTKVIEWWINEXPORT __declspec(dllexport)
  #else
    #define PYVINECOPULIB_QTVTKVIEWWINEXPORT __declspec(dllimport)
  #endif
#else
  #define PYVINECOPULIB_QTVTKVIEWWINEXPORT
#endif

#endif
