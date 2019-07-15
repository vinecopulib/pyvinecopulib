/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvQtVTKModelWin32ExportHeader_h
#define pvQtVTKModelWin32ExportHeader_h

/**
* \file pvQtVTKModelWin32ExportHeader.h
* \brief Header to sort Windows dllexport/dllimport.
*/

#if (defined(_WIN32) || defined(WIN32)) && !defined(PYVINECOPULIB_STATIC)
  #ifdef PYVINECOPULIB_QTVTKMODEL_WINDOWS_EXPORT
    #define PYVINECOPULIB_QTVTKMODELWINEXPORT __declspec(dllexport)
  #else
    #define PYVINECOPULIB_QTVTKMODELWINEXPORT __declspec(dllimport)
  #endif
#else
  #define PYVINECOPULIB_QTVTKMODELWINEXPORT
#endif

#endif
