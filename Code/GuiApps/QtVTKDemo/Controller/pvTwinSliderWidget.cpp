/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/
#include "pvTwinSliderWidget.h"

#include <cassert>

namespace pv
{

//-----------------------------------------------------------------------------
TwinSliderWidget::TwinSliderWidget(QWidget* parent)
: QWidget(parent)
{
  setupUi(this);

  bool ok = false;
  ok = connect(m_MinSpinBox, SIGNAL(valueChanged(int)), m_MinScrollBar, SLOT(setValue(int)));
  assert(ok);
  ok = connect(m_MinSpinBox, SIGNAL(valueChanged(int)), this, SLOT(OnLowValueChanged(int)));
  assert(ok);
  ok = connect(m_MaxSpinBox, SIGNAL(valueChanged(int)), m_MaxScrollBar, SLOT(setValue(int)));
  assert(ok);
  ok = connect(m_MaxSpinBox, SIGNAL(valueChanged(int)), this, SLOT(OnHighValueChanged(int)));
  assert(ok);
  ok = connect(m_MinScrollBar, SIGNAL(valueChanged(int)), m_MinSpinBox, SLOT(setValue(int)));
  assert(ok);
  ok = connect(m_MinScrollBar, SIGNAL(valueChanged(int)), this, SLOT(OnLowValueChanged(int)));
  assert(ok);
  ok = connect(m_MaxScrollBar, SIGNAL(valueChanged(int)), m_MaxSpinBox, SLOT(setValue(int)));
  assert(ok);
  ok = connect(m_MaxScrollBar, SIGNAL(valueChanged(int)), this, SLOT(OnHighValueChanged(int)));
  assert(ok);
}


//-----------------------------------------------------------------------------
TwinSliderWidget::~TwinSliderWidget()
{
}


//-----------------------------------------------------------------------------
void TwinSliderWidget::SetValues(int lowValue, int highValue)
{
  int low = lowValue;
  int high = highValue;

  if (low > high)
  {
    int tmp = low;
    low = high;
    high = tmp;
  }

  if (low < this->GetMin())
  {
    this->SetRange(low, this->GetMax());
  }
  if (high > this->GetMax())
  {
    this->SetRange(this->GetMin(), high);
  }

  m_MinSpinBox->setValue(low);
  m_MaxSpinBox->setValue(high);
}


//-----------------------------------------------------------------------------
void TwinSliderWidget::SetRange(int minValue, int maxValue)
{
  int min = minValue;
  int max = maxValue;

  if (min > max)
  {
    int tmp = min;
    min = max;
    max = tmp;
  }

  m_MinSpinBox->setMinimum(min);
  m_MinSpinBox->setMaximum(max);
  m_MaxSpinBox->setMinimum(min);
  m_MaxSpinBox->setMaximum(max);
  m_MinScrollBar->setMinimum(min);
  m_MinScrollBar->setMaximum(max);
  m_MaxScrollBar->setMinimum(min);
  m_MaxScrollBar->setMaximum(max);
}


//-----------------------------------------------------------------------------
int TwinSliderWidget::GetMin() const
{
  return m_MinSpinBox->minimum();
}


//-----------------------------------------------------------------------------
int TwinSliderWidget::GetMax() const
{
  return m_MaxSpinBox->maximum();
}


//-----------------------------------------------------------------------------
int TwinSliderWidget::GetLow() const
{
  return m_MinSpinBox->value();
}


//-----------------------------------------------------------------------------
int TwinSliderWidget::GetHigh() const
{
  return m_MaxSpinBox->value();
}


//-----------------------------------------------------------------------------
void TwinSliderWidget::OnLowValueChanged(int low)
{
  if (low > this->GetHigh())
  {
    m_MaxSpinBox->setValue(low+1);
    m_MaxScrollBar->setValue(low+1);
  }
  emit ValuesChanged(m_MinSpinBox->value(), m_MaxSpinBox->value());
}


//-----------------------------------------------------------------------------
void TwinSliderWidget::OnHighValueChanged(int high)
{
  if (high < this->GetLow())
  {
    m_MinSpinBox->setValue(high-1);
    m_MinScrollBar->setValue(high-1);
  }
  emit ValuesChanged(m_MinSpinBox->value(), m_MaxSpinBox->value());
}

} // end namespace
