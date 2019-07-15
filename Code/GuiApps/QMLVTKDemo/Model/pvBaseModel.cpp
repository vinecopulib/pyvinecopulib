/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "pvBaseModel.h"

#include <QQuickWindow>

namespace pv
{

//-----------------------------------------------------------------------------
BaseModel::BaseModel()
{
  connect(this, &QQuickItem::windowChanged, this, &BaseModel::handleWindowChanged);
}


//-----------------------------------------------------------------------------
void BaseModel::handleWindowChanged(QQuickWindow *win)
{
  if(window())
  {
    disconnect(window(), 0, this, 0);
  }

  if (win)
  {
    connect(win, &QQuickWindow::beforeSynchronizing, this, &BaseModel::sync, Qt::DirectConnection);
    connect(win, &QQuickWindow::sceneGraphInvalidated, this, &BaseModel::cleanup, Qt::DirectConnection);

    // If we allow QML to do the clearing, they would
    // clear what we paint and nothing would show.
    win->setClearBeforeRendering(false);
  }
}


//-----------------------------------------------------------------------------
void BaseModel::cleanup()
{
  this->InternalCleanup();
}


//-----------------------------------------------------------------------------
void BaseModel::sync()
{
  this->InternalSync();
}

} // end namespace
