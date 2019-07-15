/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvMainWindow_h
#define pvMainWindow_h

#include "ui_pvMainWindow.h"
#include <pvVolumeRenderingModel.h>

#include <QMainWindow>

namespace pv
{

class VTKViewWidget;

/**
* \class MainWindow
* \brief Demo Widget provides main window, and connects it to Model.
*/
class MainWindow : public QMainWindow, public Ui_MainWindow
{
  Q_OBJECT

public:

  MainWindow(pv::VolumeRenderingModel* model);
  virtual ~MainWindow();
  void ConnectRenderer();

private slots:

  void OnFileOpen();

private:

  pv::VolumeRenderingModel* m_Model;
  pv::VTKViewWidget*        m_Viewer;

}; // end class

} // end namespace

#endif
