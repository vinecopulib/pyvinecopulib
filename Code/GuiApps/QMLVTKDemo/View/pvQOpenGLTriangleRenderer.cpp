/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "pvQOpenGLTriangleRenderer.h"
#include <QtMath>
#include <iostream>

namespace pv
{

static const char *vertexShaderSource =
    "#version 150 core\n"
    "in vec2 position;\n"
    "in vec3 color;\n"
    "out vec3 Color\n;"
    "uniform mat4 projMatrix;\n"
    "uniform mat4 mvMatrix;\n"
    "void main()\n"
    "{\n"
    "  Color = color;"
    "  gl_Position = projMatrix * mvMatrix * vec4(position, 0.0, 1.0);\n"
    "}\n";

static const char *fragmentShaderSource =
    "#version 150 core\n"
    "in vec3 Color;\n"
    "out vec4 outColor;\n"
    "void main()\n"
    "{\n"
    "  outColor = vec4(Color, 1.0);;\n"
    "}\n";

//-----------------------------------------------------------------------------
QOpenGLTriangleRenderer::QOpenGLTriangleRenderer()
: m_Erase(true)
, m_Degrees(0)
, m_Window(nullptr)
, m_Program(nullptr)
, m_TriangleData(nullptr)
, m_TriangleDataDirty(false)
{
}


//-----------------------------------------------------------------------------
QOpenGLTriangleRenderer::~QOpenGLTriangleRenderer()
{
  if (m_Program != nullptr)
  {
    m_Program->bind();
    m_VBO.destroy();
    m_Program->release();
    delete m_Program;
    m_Program = nullptr;
  }
}


//-----------------------------------------------------------------------------
void QOpenGLTriangleRenderer::SetEraseBeforeVTKRendering(bool b)
{
  m_Erase = b;
}


//-----------------------------------------------------------------------------
void QOpenGLTriangleRenderer::setTriangleData(QVector<float>* data)
{
  m_TriangleData = data;
  m_TriangleDataDirty = true;
}


//-----------------------------------------------------------------------------
void QOpenGLTriangleRenderer::setViewportSize(const QSize &size)
{
  m_ViewportSize = size;
}


//-----------------------------------------------------------------------------
void QOpenGLTriangleRenderer::setWindow(QQuickWindow *window)
{
  m_Window = window;
}


//-----------------------------------------------------------------------------
void QOpenGLTriangleRenderer::paint()
{
  initializeOpenGLFunctions();

  if (!m_Program)
  {
    m_Program = new QOpenGLShaderProgram();
    m_Program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
    m_Program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);
    m_Program->bindAttributeLocation("position", 0);
    m_Program->bindAttributeLocation("color", 1);
    m_Program->link();
    m_Program->bind();
    m_ProjMatrixLoc = m_Program->uniformLocation("projMatrix");
    m_ModelViewMatrixLoc = m_Program->uniformLocation("mvMatrix");
    m_Program->release();
  }

  if (m_TriangleDataDirty)
  {
    if (!m_VAO.isCreated())
    {
      m_VAO.create();
    }

    QOpenGLVertexArrayObject::Binder vaoBinder(&m_VAO);
    if (!m_VBO.isCreated())
    {
      m_VBO.create();
      m_VBO.bind();
      m_VBO.allocate(m_TriangleData->data(), m_TriangleData->size() * sizeof(float));
      m_VBO.release();
    }

    m_Program->bind();
    m_VBO.bind();

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));

    m_VBO.release();
    m_Program->release();

    m_TriangleDataDirty = false;
  }

  glViewport(0, 0, m_ViewportSize.width(), m_ViewportSize.height());
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);

  if (m_Erase)
  {
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  }

  m_ModelViewMatrix.setToIdentity();
  m_ModelViewMatrix.rotate(m_Degrees, 0, 0, 1);

  m_CameraMatrix.setToIdentity();
  m_CameraMatrix.translate(0, 0, -3);

  m_ProjMatrix.setToIdentity();
  m_ProjMatrix.perspective(80.0f, GLfloat(m_ViewportSize.width()) / m_ViewportSize.height(), 0.01f, 100.0f);

  QOpenGLVertexArrayObject::Binder vaoBinder(&m_VAO);
  m_Program->bind();
  m_Program->setUniformValue(m_ProjMatrixLoc, m_ProjMatrix);
  m_Program->setUniformValue(m_ModelViewMatrixLoc, m_CameraMatrix * m_ModelViewMatrix);

  glDrawArrays(GL_TRIANGLES, 0, 3);

  m_Program->release();
  m_Window->resetOpenGLState();
}

} // end namespace
