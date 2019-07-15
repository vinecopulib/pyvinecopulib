/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "pvCubeModel.h"
#include <pvQQuickVTKView.h>
#include <vtkProperty.h>

namespace pv
{

//-----------------------------------------------------------------------------
CubeModel::CubeModel()
: m_Connected(false)
, m_Degrees(0)
, m_CubeSource(nullptr)
, m_CubeMapper(nullptr)
, m_CubeActor(nullptr)
, m_Renderer(nullptr)
{
  m_CubeSource = vtkCubeSource::New();
  m_CubeSource->SetCenter(0, 0, 0);
  m_CubeSource->SetXLength(1);
  m_CubeSource->SetYLength(1);
  m_CubeSource->SetZLength(1);
  m_CubeMapper = vtkPolyDataMapper::New();
  m_CubeMapper->SetInputConnection(m_CubeSource->GetOutputPort());
  m_CubeActor = vtkActor::New();
  m_CubeActor->GetProperty()->SetColor(1, 0, 0);
  m_CubeActor->GetProperty()->SetOpacity(0.5);
  m_CubeActor->SetMapper(m_CubeMapper);
  m_Renderer = vtkRenderer::New();
  m_Renderer->AddActor(m_CubeActor);
}


//-----------------------------------------------------------------------------
void CubeModel::setDegrees(qreal d)
{
  if (d == m_Degrees)
  {
    return;
  }
  m_Degrees = d;
  m_CubeActor->RotateY(m_Degrees); // this will be cumulative, but it just needs to move.
  if (window())
  {
    window()->update();
  }
}


//-----------------------------------------------------------------------------
void CubeModel::InternalCleanup()
{
  QQuickVTKView* w = dynamic_cast<QQuickVTKView*>(window());
  w->RemoveRenderer(m_Renderer);
}


//-----------------------------------------------------------------------------
void CubeModel::InternalSync()
{
  if (!m_Connected)
  {
    QQuickVTKView* w = dynamic_cast<QQuickVTKView*>(window());
    w->AddRenderer(m_Renderer);
    m_Connected = true;
  }
}


} // end namespace
