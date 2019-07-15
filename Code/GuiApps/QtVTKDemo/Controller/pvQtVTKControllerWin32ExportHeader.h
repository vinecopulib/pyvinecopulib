/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvQtVTKControllerWin32ExportHeader_h
#define pvQtVTKControllerWin32ExportHeader_h

/**
* \file pvQtVTKControllerWin32ExportHeader.h
* \brief Header to sort Windows dllexport/dllimport.
*/

#if (defined(_WIN32) || defined(WIN32)) && !defined(PYVINECOPULIB_STATIC)
  #ifdef PYVINECOPULIB_QTVTKCONTROLLER_WINDOWS_EXPORT
    #define PYVINECOPULIB_QTVTKCONTROLLERWINEXPORT __declspec(dllexport)
  #else
    #define PYVINECOPULIB_QTVTKCONTROLLERWINEXPORT __declspec(dllimport)
  #endif
#else
  #define PYVINECOPULIB_QTVTKCONTROLLERWINEXPORT
#endif

#endif
