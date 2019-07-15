/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "pvOpenGLWidget.h"
#include <QCoreApplication>
#include <QDebug>

namespace pv
{

bool OpenGLWidget::m_IsTransparent = false;

static const char *vertexShaderSource =
    "#version 150 core\n"
    "in vec2 position;\n"
    "in vec3 color;\n"
    "out vec3 Color\n;"
    "void main()\n"
    "{\n"
    "  Color = color;"
    "  gl_Position = vec4(position, 0.0, 1.0);\n"
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
OpenGLWidget::OpenGLWidget(QWidget *parent)
: QOpenGLWidget(parent)
, m_IsCore(false)
{
  m_IsCore = QSurfaceFormat::defaultFormat().profile() == QSurfaceFormat::CoreProfile;
  if (m_IsTransparent)
  {
    QSurfaceFormat fmt = format();
    fmt.setAlphaBufferSize(8);
    setFormat(fmt);
  }
}


//-----------------------------------------------------------------------------
OpenGLWidget::~OpenGLWidget()
{
  cleanup();
}


//-----------------------------------------------------------------------------
QSize OpenGLWidget::minimumSizeHint() const
{
  return QSize(50, 50);
}


//-----------------------------------------------------------------------------
QSize OpenGLWidget::sizeHint() const
{
  return QSize(400, 400);
}


//-----------------------------------------------------------------------------
void OpenGLWidget::cleanup()
{
  this->makeCurrent();

  glDeleteProgram(m_ShaderProgram);
  glDeleteShader(m_FragmentShader);
  glDeleteShader(m_VertexShader);
  glDeleteBuffers(1, &m_VBO);
  glDeleteVertexArrays(1, &m_VAO);

  this->doneCurrent();
}


//-----------------------------------------------------------------------------
void OpenGLWidget::initializeGL()
{
  // In this example the widget's corresponding top-level window can change
  // several times during the widget's lifetime. Whenever this happens, the
  // QOpenGLWidget's associated context is destroyed and a new one is created.
  // Therefore we have to be prepared to clean up the resources on the
  // aboutToBeDestroyed() signal, instead of the destructor. The emission of
  // the signal will be followed by an invocation of initializeGL() where we
  // can recreate all resources.
  connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &OpenGLWidget::cleanup);

  initializeOpenGLFunctions();

  glGenVertexArrays(1, &m_VAO);
  glBindVertexArray(m_VAO);
  glGenBuffers(1, &m_VBO);

  GLfloat vertices[] = {
     0.0f,  0.5f, 1.0f, 0.0f, 0.0f, // Vertex 1: Red
     0.5f, -0.5f, 0.0f, 1.0f, 0.0f, // Vertex 2: Green
    -0.5f, -0.5f, 0.0f, 0.0f, 1.0f  // Vertex 3: Blue
  };

  glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

  m_VertexShader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(m_VertexShader, 1, &vertexShaderSource, NULL);
  glCompileShader(m_VertexShader);

  m_FragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(m_FragmentShader, 1, &fragmentShaderSource, NULL);
  glCompileShader(m_FragmentShader);

  m_ShaderProgram = glCreateProgram();
  glAttachShader(m_ShaderProgram, m_VertexShader);
  glAttachShader(m_ShaderProgram, m_FragmentShader);
  glBindFragDataLocation(m_ShaderProgram, 0, "outColor");
  glLinkProgram(m_ShaderProgram);
  glUseProgram(m_ShaderProgram);

  m_PositionAttribute = glGetAttribLocation(m_ShaderProgram, "position");
  glEnableVertexAttribArray(m_PositionAttribute);
  glVertexAttribPointer(m_PositionAttribute, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), 0);
  m_ColourAttribute = glGetAttribLocation(m_ShaderProgram, "color");
  glEnableVertexAttribArray(m_ColourAttribute);
  glVertexAttribPointer(m_ColourAttribute, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)(2*sizeof(float)));
}


//-----------------------------------------------------------------------------
void OpenGLWidget::paintGL()
{
  const qreal retinaScale = this->devicePixelRatio();
  glViewport(0, 0, width() * retinaScale, height() * retinaScale);

  glClearColor(0.0f, 0.0f, 0.0f, m_IsTransparent ? 0 : 1);
  glClear(GL_COLOR_BUFFER_BIT);

  glDrawArrays(GL_TRIANGLES, 0, 3);
}


//-----------------------------------------------------------------------------
void OpenGLWidget::resizeGL(int w, int h)
{
}


} // end namespace
