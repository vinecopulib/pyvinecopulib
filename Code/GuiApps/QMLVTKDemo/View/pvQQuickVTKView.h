/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvQQuickVTKView_h
#define pvQQuickVTKView_h

#include <QQuickView>
#include <QMutex>
#include <vtkSmartPointer.h>
#include <vtkExternalOpenGLRenderWindow.h>
#include <vtkInteractorStyleMultiTouchCamera.h>

class vtkRenderer;
class vtkEventQtSlotConnect;
class QVTKInteractorAdapter;
class QVTKInteractor;

namespace pv
{

/**
 * \class QQuickVTKView
 * \brief Renders a VTK scene into a QQuickView
 */
class QQuickVTKView : public QQuickView {

  Q_OBJECT

public:

  virtual ~QQuickVTKView();
  QQuickVTKView(QWindow * parent = 0);
  QQuickVTKView(QQmlEngine * engine, QWindow * parent);
  QQuickVTKView(const QUrl & source, QWindow * parent = 0);

  void AddRenderer(vtkRenderer* r);
  void RemoveRenderer(vtkRenderer *r);
  void SetEraseBeforeVTKRendering(bool b);
  void SetEnabled(bool isEnabled);

protected:

  virtual bool event(QEvent* e);
  virtual void mousePressEvent(QMouseEvent* event);
  virtual void mouseMoveEvent(QMouseEvent* event);
  virtual void mouseReleaseEvent(QMouseEvent* event);
  virtual void mouseDoubleClickEvent(QMouseEvent* event);
  virtual void keyPressEvent(QKeyEvent* event);
  virtual void keyReleaseEvent(QKeyEvent* event);
  virtual void enterEvent(QEvent*);
  virtual void leaveEvent(QEvent*);
  virtual void wheelEvent(QWheelEvent*);

signals:

  /**
   * \brief Triggered before VTK renders.
   * Must be connected using Qt::DirectConnection.
   */
  void beforeVTKRendering();

  /**
   * \brief Triggered after VTK renders.
   * Must be connected using Qt::DirectConnection.
   */
  void afterVTKRendering();

private slots:

  void Render();

private:

  bool ProcessEvent(QEvent* e);

  QMutex                                              m_Mutex;
  vtkSmartPointer<vtkExternalOpenGLRenderWindow>      m_VTKRenderWindow;
  vtkSmartPointer<QVTKInteractor>                     m_VTKRenderWindowInteractor;
  vtkSmartPointer<vtkInteractorStyleMultiTouchCamera> m_VTKInteractorStyleMultiTouchCamera;
  vtkSmartPointer<vtkEventQtSlotConnect>              m_EventSlotConnector;
  QVTKInteractorAdapter*                              m_VTKInteractorAdapter;
  bool                                                m_EraseBeforeVTKRendering;
  void Init();

};

} // end namespace

#endif
