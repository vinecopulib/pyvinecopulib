/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvVTKViewWidget_h
#define pvVTKViewWidget_h

#include "pvQtVTKViewWin32ExportHeader.h"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>

#ifdef BUILD_VTK_OpenGL
#include <QVTKWidget2.h>
#else
#include <QVTKOpenGLWidget.h>
#endif

namespace pv
{

/**
* \class VTKViewWidget
* \brief Demo Widget to provide a standard VTK window within a Qt widget.
*
* Note, as of issue #16, we show how to use either VTKs OpenGL backend
* or VTKs OpenGL2 backend, which are fundamentally different, and their
* usage is also linked to which version of Qt you are using. In practice
* remember that this is an example project. You would NOT want to
* support all these options. Pick the most up to date version you can.
*
* On June 20th 2014, VTK merged changes to their master such that you could
* choose a backend:
*
* \li OpenGL - uses the traditional OpenGL 1.1 fixed pipeline calls.
* \li OpenGL2 - assumes a minimum of OpenGL API version of 2.1 and uses shaders.
*
* So, as of VTK version 6.3.0 you have the choice of either. However, the vtkRenderWindow
* must be integrated within a correct choice of Qt window which complicated matters somewhat.
*
* Within VTK, there have been different widgets to integrate a
* vtkRenderWindow with a QWidget.
*
* \li QVTKWidget, derived from QWidget, uses QPaintEngine to draw on screen.
* \li QVTKWidget2, derived from QGLWidget, and directly displays the OpenGL rendering that is contained within vtkRenderWindow, but QGLWidget is deprecated as of Qt 5.4.
* \li QVTKOpenGLWidget, derived from QOpenGLWidget, added on Jan 27th 2017 after having been developed in ParaView, available from VTK 7.1.1, and the docs in the header file says it targets Qt 5.5 and above.
*
* This project only supports Qt5 and above, so we prefer 2 choices:
*
* \li If Qt >= 5.5.0, then use OpenGL2 backend, which means QVTKOpenGLWidget
* \li If Qt < 5.5.0, then use OpenGL backend, which means QVTKWidget or QVTKWidget2, but lets chose QVTKWidget2.
*
* Also, if you look in the header file for QVTKOpenGLWidget in the VTK source code, you will seed
* more instructions for use. These have been implemented in our main.cpp.
*/
class PYVINECOPULIB_QTVTKVIEWWINEXPORT VTKViewWidget
#ifdef BUILD_VTK_OpenGL
    : public QVTKWidget2
#else
    : public QVTKOpenGLWidget
#endif
{
  Q_OBJECT

public:

  VTKViewWidget(QWidget* parent);
  virtual ~VTKViewWidget();

  void AddRenderer(vtkRenderer* r);

public slots:

  void Render();

private:

#ifdef BUILD_VTK_OpenGL2
  vtkSmartPointer<vtkGenericOpenGLRenderWindow> m_RenderWindow;
#endif

}; // end class

} // end namespace

#endif
