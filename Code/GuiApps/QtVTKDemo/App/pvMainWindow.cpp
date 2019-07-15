/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/
#include "pvMainWindow.h"
#include <pvExceptionMacro.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <QFileDialog>
#include <QStandardPaths>

#include <cassert>

namespace pv
{

//-----------------------------------------------------------------------------
MainWindow::MainWindow(pv::VolumeRenderingModel* model)
{
  if (model == nullptr)
  {
    pvExceptionThrow() << "Model is null.";
  }

  setupUi(this);
  setCentralWidget(m_CentralWidget);
  this->setWindowIcon(QIcon(":/QtVTK/icon.png"));

  m_Model = model;

  bool ok = false;

  // Connect main window to model.
  ok = connect(actionOpen, SIGNAL(triggered()), this, SLOT(OnFileOpen()));
  assert(ok);

  // Connect widgets to model
  ok = connect(m_CentralWidget, SIGNAL(WindowValuesChanged(int, int)), m_Model, SLOT(SetIntensityWindow(int, int)));
  assert(ok);
  ok = connect(m_CentralWidget, SIGNAL(DoSomethingPressed()), m_Model, SLOT(DoSomethingPressed()));
  assert(ok);

  // Connect model to widgets.
  ok = connect(m_Model, SIGNAL(ImageLoaded(int, int)), m_CentralWidget, SLOT(SetIntensityRange(int, int)));
  assert(ok);
  ok = connect(m_Model, SIGNAL(Modified()), m_CentralWidget->GetVTKViewWidget(), SLOT(Render()));
  assert(ok);
}


//-----------------------------------------------------------------------------
MainWindow::~MainWindow()
{
}


//-----------------------------------------------------------------------------
void MainWindow::ConnectRenderer()
{
  m_Viewer = m_CentralWidget->GetVTKViewWidget();
  m_Viewer->AddRenderer(m_Model->GetRenderer());
}


//-----------------------------------------------------------------------------
void MainWindow::OnFileOpen()
{

  QStringList paths;
  paths = QStandardPaths::standardLocations(QStandardPaths::HomeLocation);
  assert(paths.size() == 1);
  QString path = paths[0];

  QString dirName = QFileDialog::getExistingDirectory(
        this, "Choose a directory", path,
        QFileDialog::ShowDirsOnly);

  if (!dirName.isEmpty())
  {
    m_Model->LoadDirectory(dirName.toStdString());
  }
}

} // end namespace
