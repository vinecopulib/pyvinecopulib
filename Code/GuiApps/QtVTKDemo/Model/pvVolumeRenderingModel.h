/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef pvVolumeRenderingModel_h
#define pvVolumeRenderingModel_h

#include "pvQtVTKModelWin32ExportHeader.h"

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkDICOMImageReader.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkVolume.h>
#include <QObject>

namespace pv
{

/**
* \class VolumeRenderingModel
* \brief Demo Model, to contain VTK pipelines.
*
* Intended to demonstrate that this class knows nothing about the View.
* To further illustrate this point, this class is not a QWidget,
* it derives from QObject. This class is controlled via Qt signals and slots.
*
* It would be possible to make some of the rendering pipeline part of the View.
* However, the VTK rendering pipeline contains parameters, which must be set,
* and these parameters form part of the data-model. So, for simplicity's sake
* I have kept it all in the model. The View layer therefore just renders
* the actors, so should remain quite thin.
*/
  class PYVINECOPULIB_QTVTKMODELWINEXPORT VolumeRenderingModel : public QObject
{
  Q_OBJECT

public:

  VolumeRenderingModel();
  virtual ~VolumeRenderingModel();

  /**
   * \brief Load DICOM directory into the pipeline.
   */
  void LoadDirectory(const std::string& dirName);

  /**
   * \brief Get hold of the renderer.
   */
  vtkRenderer* GetRenderer() const;

signals:

  void ImageLoaded(int minValue, int maxValue);
  void Modified();

public slots:

  void SetIntensityWindow(int minValue, int maxValue);
  void DoSomethingPressed();

private:

  vtkSmartPointer<vtkDICOMImageReader>      m_ImageReader;
  vtkSmartPointer<vtkSmartVolumeMapper>     m_VolumeMapper;
  vtkSmartPointer<vtkColorTransferFunction> m_ColorTransferFunction;
  vtkSmartPointer<vtkPiecewiseFunction>     m_OpacityTransferFunction;
  vtkSmartPointer<vtkVolumeProperty>        m_VolumeProperty;
  vtkSmartPointer<vtkVolume>                m_Volume;
  vtkSmartPointer<vtkRenderer>              m_Renderer;

}; // end class

} // end namespace

#endif
