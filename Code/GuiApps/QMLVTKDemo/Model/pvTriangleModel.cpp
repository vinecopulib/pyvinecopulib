/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "pvTriangleModel.h"

#include <pvQOpenGLTriangleRenderer.h>
#ifdef BUILD_VTK // for demo purposes only. In practice I wouldnt mix VTK with other OpenGL.
#include <pvQQuickVTKView.h>
#endif

namespace pv
{

//-----------------------------------------------------------------------------
TriangleModel::TriangleModel()
: m_Degrees(0)
, m_Renderer(nullptr)
, m_TriangleData({
                 0.0f,  0.5f, 1.0f, 0.0f, 0.0f, // Vertex 1: position, red
                 0.5f, -0.5f, 0.0f, 1.0f, 0.0f, // Vertex 2: position, green
                -0.5f, -0.5f, 0.0f, 0.0f, 1.0f  // Vertex 3: position, blue
                 })
{
}


//-----------------------------------------------------------------------------
void TriangleModel::setDegrees(qreal d)
{
  if (d == m_Degrees)
  {
    return;
  }
  m_Degrees = d;
  if (window())
  {
    window()->update();
  }
}


//-----------------------------------------------------------------------------
void TriangleModel::InternalCleanup()
{
  if (m_Renderer)
  {
    delete m_Renderer;
    m_Renderer = nullptr;
  }
}


//-----------------------------------------------------------------------------
void TriangleModel::InternalSync()
{
  if (!m_Renderer)
  {
    m_Renderer = new QOpenGLTriangleRenderer();
    m_Renderer->setTriangleData(&m_TriangleData);
    m_Renderer->setWindow(window());
    m_Renderer->SetEraseBeforeVTKRendering(false);
    connect(dynamic_cast<QQuickVTKView*>(window()), &QQuickVTKView::afterVTKRendering,
            m_Renderer, &QOpenGLTriangleRenderer::paint, Qt::DirectConnection);
  }
  m_Renderer->setViewportSize(window()->size() * window()->devicePixelRatio());
  m_Renderer->setDegrees(m_Degrees);
}

} // end namespace
