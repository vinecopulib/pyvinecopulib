/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvOpenGLWidget_h
#define pvOpenGLWidget_h

#include <QOpenGLWidget>
#include <QOpenGLFunctions>

namespace pv
{

/**
 * \class OpenGLWidget
 * \brief Demo Widget to show how to setup an OpenGL window, derived from Qt's QOpenGLWidget.
 */
class OpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
  Q_OBJECT

public:

  OpenGLWidget(QWidget *parent = 0);
  ~OpenGLWidget();

  static bool isTransparent() { return m_IsTransparent; }
  static void setTransparent(bool t) { m_IsTransparent = t; }

  QSize minimumSizeHint() const override;
  QSize sizeHint() const override;

public slots:

  void cleanup();

protected:

  void initializeGL() override;
  void paintGL() override;
  void resizeGL(int width, int height) override;

private:

  static bool m_IsTransparent;
  bool        m_IsCore;

  GLuint      m_VAO;
  GLuint      m_VBO;
  GLuint      m_VertexShader;
  GLuint      m_FragmentShader;
  GLuint      m_ShaderProgram;
  GLuint      m_PositionAttribute;
  GLuint      m_ColourAttribute;
};

} // end namespace

#endif
