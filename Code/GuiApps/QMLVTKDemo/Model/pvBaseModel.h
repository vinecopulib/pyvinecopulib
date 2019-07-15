/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvBaseModel_h
#define pvBaseModel_h

#include <QtQuick/QQuickItem>

namespace pv
{

/**
 * \class BaseModel
 * \brief Base class for our model classes, each derived from QQuickItem.
 */
class BaseModel : public QQuickItem
{
  Q_OBJECT

public:

  BaseModel();

public slots:

  /**
   * \name Render Thread Methods
   */
  ///@{

  /// \brief Called when the current window emits QQuickWindow::beforeSynchronizing.
  void sync();

  /// \brief Called when the current window emits QQuickWindow::sceneGraphInvalidated.
  void cleanup();

  ///@}

protected:

  virtual void InternalCleanup() {}
  virtual void InternalSync() {}

private slots:

  void handleWindowChanged(QQuickWindow *win);

};

} // end namespace

#endif
