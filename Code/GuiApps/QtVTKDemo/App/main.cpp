/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include <QVTKApplication.h>
#include "PyVinecopulibConfigure.h"
#include "pvMainWindow.h"
#include <pvVolumeRenderingModel.h>
#include <QScopedPointer>

#ifdef BUILD_VTK_OpenGL2
#include <QSurfaceFormat>
#include <QVTKOpenGLWidget.h>
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkNew.h"
#endif

int main(int argc, char** argv)
{

#ifdef BUILD_VTK_OpenGL2
  vtkOpenGLRenderWindow::SetGlobalMaximumNumberOfMultiSamples(0);
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
#endif

  QVTKApplication app(argc, argv);
  app.setOrganizationName("UCL");
  app.setApplicationName("QtVTKDemo");
  app.setApplicationVersion(QString(PYVINECOPULIB_VERSION_STRING));

  QScopedPointer<pv::VolumeRenderingModel> mb(new pv::VolumeRenderingModel());

  pv::MainWindow mainWin(mb.data());

  mainWin.show();
  mainWin.ConnectRenderer();
  mainWin.showMaximized();

  return app.exec();
}
