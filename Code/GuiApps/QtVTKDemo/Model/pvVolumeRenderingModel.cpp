/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/
#include "pvVolumeRenderingModel.h"

#include <vtkImageData.h>
#include <vtkImageHistogramStatistics.h>

#include <stdexcept>
#include <cassert>

namespace pv
{

//-----------------------------------------------------------------------------
VolumeRenderingModel::VolumeRenderingModel()
{
  m_ImageReader = vtkSmartPointer<vtkDICOMImageReader>::New();

  m_VolumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
  m_VolumeMapper->SetInputConnection(m_ImageReader->GetOutputPort());

  m_ColorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
  m_OpacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();

  m_VolumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
  m_VolumeProperty->SetIndependentComponents(true);
  m_VolumeProperty->SetColor(m_ColorTransferFunction);
  m_VolumeProperty->SetScalarOpacity(m_OpacityTransferFunction);
  m_VolumeProperty->SetInterpolationTypeToLinear();

  m_Volume = vtkSmartPointer<vtkVolume>::New();
  m_Volume->SetProperty(m_VolumeProperty);
  m_Volume->SetMapper(m_VolumeMapper);

  m_Renderer = vtkSmartPointer<vtkRenderer>::New();
}


//-----------------------------------------------------------------------------
VolumeRenderingModel::~VolumeRenderingModel()
{
}


//-----------------------------------------------------------------------------
vtkRenderer* VolumeRenderingModel::GetRenderer() const
{
  return m_Renderer.GetPointer();
}


//-----------------------------------------------------------------------------
void VolumeRenderingModel::LoadDirectory(const std::string& dirName)
{
  m_ImageReader->SetDirectoryName(dirName.c_str());
  m_ImageReader->Update();

  vtkImageData *input = m_ImageReader->GetOutput();
  int *dims = input->GetDimensions();

  if (dims[0] > 1 && dims[1] > 1 && dims[2] > 1)
  {
    vtkSmartPointer<vtkImageHistogramStatistics> ihs =
      vtkSmartPointer<vtkImageHistogramStatistics>::New();
    ihs->SetInputConnection(m_ImageReader->GetOutputPort());
    ihs->Update();

    double min = ihs->GetMinimum();
    double max = ihs->GetMaximum();

    std::cout << "Loaded image:" << dirName << ", size="
              << dims[0] << ", "
              << dims[1] << ", "
              << dims[2] << ", "
              << ", range=" << min << ", " << max
              << std::endl;

    m_Renderer->AddVolume(m_Volume);
    m_Renderer->ResetCamera();

    m_ColorTransferFunction->AddRGBSegment(min, 1.0, 1.0, 1.0, max, 1.0, 1.0, 1.0 );
    this->SetIntensityWindow(min, max);

    emit ImageLoaded(min, max);
  }
}


//-----------------------------------------------------------------------------
void VolumeRenderingModel::SetIntensityWindow(int minValue, int maxValue)
{
  m_OpacityTransferFunction->AddSegment(minValue, 0.0,
                                        maxValue, 1.0 );
  m_VolumeMapper->SetBlendModeToComposite();
  m_VolumeProperty->ShadeOn();

  emit Modified();
}


//-----------------------------------------------------------------------------
void VolumeRenderingModel::DoSomethingPressed()
{
  std::cout << "DoSomethingPressed()" << std::endl;
}

} // end namespace
