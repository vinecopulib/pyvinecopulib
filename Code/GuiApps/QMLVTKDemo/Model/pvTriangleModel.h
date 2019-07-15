/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvTriangleModel_h
#define pvTriangleModel_h

#include <QVector>
#include "pvBaseModel.h"
#include <pvQOpenGLTriangleRenderer.h>

namespace pv
{

/**
 * \class TriangleModel
 * \brief Demo Model class containing data for a single, rotateable, triangle.
 */
class TriangleModel : public BaseModel
{
  Q_OBJECT
  Q_PROPERTY(qreal degrees READ degrees WRITE setDegrees)

public:

  TriangleModel();

  /**
   * \name Data Setter Methods
   */
  ///@{

  qreal degrees() const { return m_Degrees; }
  void setDegrees(qreal d);

  ///@}

protected:

  virtual void InternalCleanup();
  virtual void InternalSync();

private:

  qreal                    m_Degrees;
  QOpenGLTriangleRenderer *m_Renderer;
  QVector<float>           m_TriangleData;

};

} // end namespace

#endif
