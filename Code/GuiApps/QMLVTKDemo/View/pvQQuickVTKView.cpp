/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "pvQQuickVTKView.h"
#include <vtkRendererCollection.h>
#include <vtkEventQtSlotConnect.h>
#include <QVTKInteractorAdapter.h>
#include <QVTKInteractor.h>
#include <QMutexLocker>

namespace pv {

//-----------------------------------------------------------------------------
QQuickVTKView::~QQuickVTKView()
{
  if (m_VTKRenderWindow)
  {
    m_VTKRenderWindowInteractor->Disable();
    m_VTKRenderWindow->SetInteractor(nullptr);
#ifdef BUILD_VTK_OpenGL2
    m_VTKRenderWindow->SetReadyForRendering(false);
#endif
    m_VTKRenderWindow->RemoveAllObservers();
    vtkRendererCollection *c = m_VTKRenderWindow->GetRenderers();
    vtkRenderer *v = nullptr;
    while(v = c->GetNextItem(), v != nullptr)
    {
      m_VTKRenderWindow->RemoveRenderer(v);
    }
  }
}


//-----------------------------------------------------------------------------
QQuickVTKView::QQuickVTKView(QWindow * parent)
: QQuickView(parent)
, m_Mutex(QMutex::NonRecursive)
, m_VTKRenderWindow(nullptr)
, m_VTKRenderWindowInteractor(nullptr)
, m_VTKInteractorStyleMultiTouchCamera(nullptr)
, m_EraseBeforeVTKRendering(true)
{
  this->Init();
}


//-----------------------------------------------------------------------------
QQuickVTKView::QQuickVTKView(QQmlEngine * engine, QWindow * parent)
: QQuickView(engine, parent)
, m_Mutex(QMutex::NonRecursive)
, m_VTKRenderWindow(nullptr)
, m_VTKRenderWindowInteractor(nullptr)
, m_VTKInteractorStyleMultiTouchCamera(nullptr)
, m_EraseBeforeVTKRendering(true)
{
  this->Init();
}


//-----------------------------------------------------------------------------
QQuickVTKView::QQuickVTKView(const QUrl & source, QWindow * parent)
: QQuickView(source, parent)
, m_Mutex(QMutex::NonRecursive)
, m_VTKRenderWindow(nullptr)
, m_VTKRenderWindowInteractor(nullptr)
, m_VTKInteractorStyleMultiTouchCamera(nullptr)
, m_EraseBeforeVTKRendering(true)
{
  this->Init();
}


//-----------------------------------------------------------------------------
void QQuickVTKView::Init()
{
  m_VTKRenderWindow = vtkExternalOpenGLRenderWindow::New();

  m_VTKInteractorStyleMultiTouchCamera = vtkInteractorStyleMultiTouchCamera::New();

  m_VTKRenderWindowInteractor = QVTKInteractor::New();
  m_VTKRenderWindowInteractor->SetRenderWindow(m_VTKRenderWindow);
  m_VTKRenderWindowInteractor->SetInteractorStyle(m_VTKInteractorStyleMultiTouchCamera);
  m_VTKRenderWindowInteractor->Initialize();

  // This is to stop interactor triggering Render(),
  // as the rendering is done in the scene graph thread,
  // but the interaction is done in the gui thread.
  m_VTKRenderWindowInteractor->SetEnableRender(false);

  m_VTKInteractorAdapter = new QVTKInteractorAdapter(this);
  m_VTKInteractorAdapter->SetDevicePixelRatio(this->devicePixelRatio());

  m_EventSlotConnector = vtkSmartPointer<vtkEventQtSlotConnect>::New();

  connect(this, SIGNAL(beforeRendering()), this, SLOT(Render()), Qt::DirectConnection);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::AddRenderer(vtkRenderer* r)
{
  m_VTKRenderWindow->AddRenderer(r);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::RemoveRenderer(vtkRenderer *r)
{
  m_VTKRenderWindow->RemoveRenderer(r);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::SetEraseBeforeVTKRendering(bool b)
{
  m_EraseBeforeVTKRendering = b;
  m_VTKRenderWindow->SetErase(m_EraseBeforeVTKRendering);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::Render()
{
  emit beforeVTKRendering();

  QMutexLocker lock(&m_Mutex);

  m_VTKRenderWindow->Render();
  emit afterVTKRendering();
}


//-----------------------------------------------------------------------------
void QQuickVTKView::SetEnabled(bool isEnabled)
{
  if (isEnabled)
  {
    m_VTKRenderWindowInteractor->Enable();
  }
  else
  {
    m_VTKRenderWindowInteractor->Disable();
  }
}


//-----------------------------------------------------------------------------
bool QQuickVTKView::event(QEvent *e)
{
  if(e->type() == QEvent::TouchBegin ||
     e->type() == QEvent::TouchUpdate ||
     e->type() == QEvent::TouchEnd
    )
  {
    QMutexLocker lock(&m_Mutex);
    if(m_VTKRenderWindow 
#ifdef BUILD_VTK_OpenGL2
      && m_VTKRenderWindow->GetReadyForRendering()
#endif
      )
    {
      m_VTKInteractorAdapter->ProcessEvent(e, m_VTKRenderWindow->GetInteractor());
      if (e->isAccepted())
      {
        return true;
      }
    }
  }
  return QQuickView::event(e);
}


//-----------------------------------------------------------------------------
bool QQuickVTKView::ProcessEvent(QEvent* e)
{
  QMutexLocker lock(&m_Mutex);
  if (m_VTKRenderWindow 
#ifdef BUILD_VTK_OpenGL2
    && m_VTKRenderWindow->GetReadyForRendering()
#endif
    )
  {
    return m_VTKInteractorAdapter->ProcessEvent(e, m_VTKRenderWindowInteractor);
  }
  return false;
}


//-----------------------------------------------------------------------------
void QQuickVTKView::mousePressEvent(QMouseEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::mouseMoveEvent(QMouseEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::mouseReleaseEvent(QMouseEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::mouseDoubleClickEvent(QMouseEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::keyPressEvent(QKeyEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::keyReleaseEvent(QKeyEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::enterEvent(QEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::leaveEvent(QEvent* event)
{
  this->ProcessEvent(event);
}


//-----------------------------------------------------------------------------
void QQuickVTKView::wheelEvent(QWheelEvent* event)
{
  this->ProcessEvent(event);
}

//-----------------------------------------------------------------------------
} // end namespace
