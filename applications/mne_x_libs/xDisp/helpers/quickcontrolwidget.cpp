//=============================================================================================================
/**
* @file     quickcontrolwidget.cpp
* @author   Lorenz Esch <Lorenz.Esch@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     June, 2015
*
* @section  LICENSE
*
* Copyright (C) 2015, Lorenz Esch and Matti Hamalainen. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that
* the following conditions are met:
*     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
*       following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
*       the following disclaimer in the documentation and/or other materials provided with the distribution.
*     * Neither the name of MNE-CPP authors nor the names of its contributors may be used
*       to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
* PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*
* @brief    Definition of the QuickControlWidget Class.
*
*/

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "quickcontrolwidget.h"


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace XDISPLIB;


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

QuickControlWidget::QuickControlWidget(QMap< qint32,float > qMapChScaling, const FiffInfo::SPtr pFiffInfo, QString name, QWidget *parent, bool bScaling, bool bProjections, bool bView, bool bFilter, bool bModalities, bool bCompensator, bool bTriggerDetection)
: RoundedEdgesWidget(parent, Qt::Window | Qt::FramelessWindowHint | Qt::WindowSystemMenuHint)
, ui(new Ui::QuickControlWidget)
, m_qMapChScaling(qMapChScaling)
, m_pFiffInfo(pFiffInfo)
, m_bScaling(bScaling)
, m_bProjections(bProjections)
, m_bView(bView)
, m_bFilter(bFilter)
, m_bModalitiy(bModalities)
, m_bCompensator(bCompensator)
, m_bTriggerDetection(bTriggerDetection)
, m_sName(name)
{
    ui->setupUi(this);

    //Init and connect hide all group (minimize) button
    ui->m_pushButton_hideAll->setText(ui->m_pushButton_hideAll->text().append(QString(" - %1").arg(m_sName)));
    connect(ui->m_pushButton_hideAll, static_cast<void (QPushButton::*)(bool)>(&QPushButton::clicked),
            this, &QuickControlWidget::toggleHideAll);

    connect(ui->m_pushButton_close, static_cast<void (QPushButton::*)(bool)>(&QPushButton::clicked),
            this, &QuickControlWidget::hide);

    //Create trigger color map
    m_qMapTriggerColor.clear();

    for(int i = 0; i<pFiffInfo->chs.size(); i++) {
        if(pFiffInfo->chs[i].kind == FIFFV_STIM_CH)
            m_qMapTriggerColor.insert(pFiffInfo->chs[i].ch_name, QColor(170,0,0));
    }

    //Create different quick control groups
    if(m_bScaling)
        createScalingGroup();
    else
        ui->m_groupBox_scaling->hide();

    if(m_bProjections)
        createProjectorGroup();
    else
        ui->m_groupBox_projections->hide();

    if(m_bView)
        createOtherGroup();
    else
        ui->m_groupBox_view->hide();

    if(m_bModalitiy)
        createModalityGroup();

    if(m_bCompensator)
        createCompensatorGroup();
    else
        ui->m_groupBox_compensators->hide();

    ui->m_groupBox_filter->hide();

    this->adjustSize();
}


//*************************************************************************************************************

QuickControlWidget::~QuickControlWidget()
{
    delete ui;
}


//*************************************************************************************************************

void QuickControlWidget::filterGroupChanged(QList<QCheckBox*> list)
{
    if(m_bFilter) {
        ui->m_groupBox_filter->show();

        m_qFilterListCheckBox.clear();

        for(int u = 0; u<list.size(); u++) {
            QCheckBox* tempCheckBox = new QCheckBox(list[u]->text());
            tempCheckBox->setChecked(list[u]->isChecked());

            connect(tempCheckBox, &QCheckBox::toggled,
                    list[u], &QCheckBox::setChecked);

            if(tempCheckBox->text() == "Activate user designed filter")
                connect(tempCheckBox, &QCheckBox::toggled,
                        this, &QuickControlWidget::userFilterToggled);

            connect(list[u], &QCheckBox::toggled,
                    tempCheckBox, &QCheckBox::setChecked);

            m_qFilterListCheckBox.append(tempCheckBox);
        }

        //Delete all widgets in the filter layout
        QGridLayout* topLayout = static_cast<QGridLayout*>(ui->m_groupBox_filter->layout());
        if(!topLayout)
           topLayout = new QGridLayout();

        QLayoutItem *child;
        while ((child = topLayout->takeAt(0)) != 0) {
            delete child->widget();
            delete child;
        }

        //Add filters
        int u = 0;

        for(u; u<m_qFilterListCheckBox.size(); u++)
            topLayout->addWidget(m_qFilterListCheckBox[u], u, 0);

        //Add push button for filter options
        m_pShowFilterOptions = new QPushButton();
//        m_pShowFilterOptions->setText("Open Filter options");
        m_pShowFilterOptions->setText("Filter options");
        m_pShowFilterOptions->setCheckable(false);
        connect(m_pShowFilterOptions, &QPushButton::clicked,
                this, &QuickControlWidget::onShowFilterOptions);

        topLayout->addWidget(m_pShowFilterOptions, u+1, 0);

        ui->m_groupBox_filter->setLayout(topLayout);

        //createViewGroup();
    } else
        ui->m_groupBox_filter->hide();
}


//*************************************************************************************************************

void QuickControlWidget::setViewParameters(double zoomFactor, int windowSize, int opactiy)
{
    ui->m_doubleSpinBox_numberVisibleChannels->setValue(zoomFactor);
    ui->m_spinBox_windowSize->setValue(windowSize);
    ui->m_horizontalSlider_opacity->setValue(opactiy);

    zoomChanged(zoomFactor);
    timeWindowChanged(windowSize);
    onOpacityChange(opactiy);
}


//*************************************************************************************************************

int QuickControlWidget::getOpacityValue()
{
    return ui->m_horizontalSlider_opacity->value();
}


//*************************************************************************************************************

int QuickControlWidget::getDistanceTimeSpacerIndex()
{
    return ui->m_comboBox_distaceTimeSpacer->currentIndex();
}


//*************************************************************************************************************

void QuickControlWidget::setDistanceTimeSpacerIndex(int index)
{
    ui->m_comboBox_distaceTimeSpacer->setCurrentIndex(index);
}


//*************************************************************************************************************

void QuickControlWidget::setNumberDetectedTriggers(int numberDetections)
{
    if(m_bTriggerDetection)
        ui->m_label_numberDetectedTriggers->setText(QString("%1").arg(numberDetections));
}


//*************************************************************************************************************

void QuickControlWidget::createScalingGroup()
{
    QGridLayout* t_pGridLayout = new QGridLayout;

    qint32 i = 0;
    //MAG
    if(m_qMapChScaling.contains(FIFF_UNIT_T))
    {
        QLabel* t_pLabelModality = new QLabel("MAG (pT)");
        t_pGridLayout->addWidget(t_pLabelModality,i,0,1,1);

        QDoubleSpinBox* t_pDoubleSpinBoxScale = new QDoubleSpinBox;
        t_pDoubleSpinBoxScale->setMinimum(0.1);
        t_pDoubleSpinBoxScale->setMaximum(50000);
        t_pDoubleSpinBoxScale->setMaximumWidth(500);
        t_pDoubleSpinBoxScale->setSingleStep(0.1);
        t_pDoubleSpinBoxScale->setDecimals(1);
        t_pDoubleSpinBoxScale->setPrefix("+/- ");
        t_pDoubleSpinBoxScale->setValue(m_qMapChScaling.value(FIFF_UNIT_T)/(1e-12));
        m_qMapScalingDoubleSpinBox.insert(FIFF_UNIT_T,t_pDoubleSpinBoxScale);
        connect(t_pDoubleSpinBoxScale,static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
                this,&QuickControlWidget::updateSpinBoxScaling);
        t_pGridLayout->addWidget(t_pDoubleSpinBoxScale,i+1,0,1,1);

        QSlider* t_pHorizontalSlider = new QSlider(Qt::Horizontal);
        t_pHorizontalSlider->setMinimum(1);
        t_pHorizontalSlider->setMaximum(5000);
        t_pHorizontalSlider->setSingleStep(1);
        t_pHorizontalSlider->setPageStep(1);
        t_pHorizontalSlider->setValue(m_qMapChScaling.value(FIFF_UNIT_T)/(1e-12)*10);
        m_qMapScalingSlider.insert(FIFF_UNIT_T,t_pHorizontalSlider);
        connect(t_pHorizontalSlider,static_cast<void (QSlider::*)(int)>(&QSlider::valueChanged),
                this,&QuickControlWidget::updateSliderScaling);
        t_pGridLayout->addWidget(t_pHorizontalSlider,i+1,1,1,1);

        i+=2;
    }

    //GRAD
    if(m_qMapChScaling.contains(FIFF_UNIT_T_M))
    {
        QLabel* t_pLabelModality = new QLabel;
        t_pLabelModality->setText("GRAD (fT/cm)");
        t_pGridLayout->addWidget(t_pLabelModality,i,0,1,1);

        QDoubleSpinBox* t_pDoubleSpinBoxScale = new QDoubleSpinBox;
        t_pDoubleSpinBoxScale->setMinimum(1);
        t_pDoubleSpinBoxScale->setMaximum(500000);
        t_pDoubleSpinBoxScale->setMaximumWidth(100);
        t_pDoubleSpinBoxScale->setSingleStep(1);
        t_pDoubleSpinBoxScale->setDecimals(1);
        t_pDoubleSpinBoxScale->setPrefix("+/- ");
        t_pDoubleSpinBoxScale->setValue(m_qMapChScaling.value(FIFF_UNIT_T_M)/(1e-15 * 100));
        m_qMapScalingDoubleSpinBox.insert(FIFF_UNIT_T_M,t_pDoubleSpinBoxScale);
        connect(t_pDoubleSpinBoxScale,static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
                this,&QuickControlWidget::updateSpinBoxScaling);
        t_pGridLayout->addWidget(t_pDoubleSpinBoxScale,i+1,0,1,1);

        QSlider* t_pHorizontalSlider = new QSlider(Qt::Horizontal);
        t_pHorizontalSlider->setMinimum(1);
        t_pHorizontalSlider->setMaximum(5000);
        t_pHorizontalSlider->setSingleStep(10);
        t_pHorizontalSlider->setPageStep(10);
        t_pHorizontalSlider->setValue(m_qMapChScaling.value(FIFF_UNIT_T_M)/(1e-15*100));
        m_qMapScalingSlider.insert(FIFF_UNIT_T_M,t_pHorizontalSlider);
        connect(t_pHorizontalSlider,static_cast<void (QSlider::*)(int)>(&QSlider::valueChanged),
                this,&QuickControlWidget::updateSliderScaling);
        t_pGridLayout->addWidget(t_pHorizontalSlider,i+1,1,1,1);

        i+=2;
    }

    //EEG
    if(m_qMapChScaling.contains(FIFFV_EEG_CH))
    {
        QLabel* t_pLabelModality = new QLabel;
        t_pLabelModality->setText("EEG (uV)");
        t_pGridLayout->addWidget(t_pLabelModality,i,0,1,1);

        QDoubleSpinBox* t_pDoubleSpinBoxScale = new QDoubleSpinBox;
        t_pDoubleSpinBoxScale->setMinimum(0.1);
        t_pDoubleSpinBoxScale->setMaximum(25000);
        t_pDoubleSpinBoxScale->setMaximumWidth(100);
        t_pDoubleSpinBoxScale->setSingleStep(0.1);
        t_pDoubleSpinBoxScale->setDecimals(1);
        t_pDoubleSpinBoxScale->setPrefix("+/- ");
        t_pDoubleSpinBoxScale->setValue(m_qMapChScaling.value(FIFFV_EEG_CH)/(1e-06));
        m_qMapScalingDoubleSpinBox.insert(FIFFV_EEG_CH,t_pDoubleSpinBoxScale);
        connect(t_pDoubleSpinBoxScale,static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
                this,&QuickControlWidget::updateSpinBoxScaling);
        t_pGridLayout->addWidget(t_pDoubleSpinBoxScale,i+1,0,1,1);

        QSlider* t_pHorizontalSlider = new QSlider(Qt::Horizontal);
        t_pHorizontalSlider->setMinimum(1);
        t_pHorizontalSlider->setMaximum(25000);
        t_pHorizontalSlider->setSingleStep(1);
        t_pHorizontalSlider->setPageStep(1);
        t_pHorizontalSlider->setValue(m_qMapChScaling.value(FIFFV_EEG_CH)/(1e-06)*10);
        m_qMapScalingSlider.insert(FIFFV_EEG_CH,t_pHorizontalSlider);
        connect(t_pHorizontalSlider,static_cast<void (QSlider::*)(int)>(&QSlider::valueChanged),
                this,&QuickControlWidget::updateSliderScaling);
        t_pGridLayout->addWidget(t_pHorizontalSlider,i+1,1,1,1);

        i+=2;
    }

    //EOG
    if(m_qMapChScaling.contains(FIFFV_EOG_CH))
    {
        QLabel* t_pLabelModality = new QLabel;
        t_pLabelModality->setText("EOG (uV)");
        t_pGridLayout->addWidget(t_pLabelModality,i,0,1,1);

        QDoubleSpinBox* t_pDoubleSpinBoxScale = new QDoubleSpinBox;
        t_pDoubleSpinBoxScale->setMinimum(0.1);
        t_pDoubleSpinBoxScale->setMaximum(102500e14);
        t_pDoubleSpinBoxScale->setMaximumWidth(100);
        t_pDoubleSpinBoxScale->setSingleStep(0.1);
        t_pDoubleSpinBoxScale->setDecimals(1);
        t_pDoubleSpinBoxScale->setPrefix("+/- ");
        t_pDoubleSpinBoxScale->setValue(m_qMapChScaling.value(FIFFV_EOG_CH)/(1e-06));
        m_qMapScalingDoubleSpinBox.insert(FIFFV_EOG_CH,t_pDoubleSpinBoxScale);
        connect(t_pDoubleSpinBoxScale,static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
                this,&QuickControlWidget::updateSpinBoxScaling);
        t_pGridLayout->addWidget(t_pDoubleSpinBoxScale,i+1,0,1,1);

        QSlider* t_pHorizontalSlider = new QSlider(Qt::Horizontal);
        t_pHorizontalSlider->setMinimum(1);
        t_pHorizontalSlider->setMaximum(25000);
        t_pHorizontalSlider->setSingleStep(1);
        t_pHorizontalSlider->setPageStep(1);
        t_pHorizontalSlider->setValue(m_qMapChScaling.value(FIFFV_EOG_CH)/(1e-06)*10);
        m_qMapScalingSlider.insert(FIFFV_EOG_CH,t_pHorizontalSlider);
        connect(t_pHorizontalSlider,static_cast<void (QSlider::*)(int)>(&QSlider::valueChanged),
                this,&QuickControlWidget::updateSliderScaling);
        t_pGridLayout->addWidget(t_pHorizontalSlider,i+1,1,1,1);

        i+=2;
    }

    //STIM
    if(m_qMapChScaling.contains(FIFFV_STIM_CH))
    {
        QLabel* t_pLabelModality = new QLabel;
        t_pLabelModality->setText("STIM");
        t_pGridLayout->addWidget(t_pLabelModality,i,0,1,1);

        QDoubleSpinBox* t_pDoubleSpinBoxScale = new QDoubleSpinBox;
        t_pDoubleSpinBoxScale->setMinimum(0.1);
        t_pDoubleSpinBoxScale->setMaximum(1000);
        t_pDoubleSpinBoxScale->setMaximumWidth(100);
        t_pDoubleSpinBoxScale->setSingleStep(0.1);
        t_pDoubleSpinBoxScale->setDecimals(1);
        t_pDoubleSpinBoxScale->setPrefix("+/- ");
        t_pDoubleSpinBoxScale->setValue(m_qMapChScaling.value(FIFFV_STIM_CH));
        m_qMapScalingDoubleSpinBox.insert(FIFFV_STIM_CH,t_pDoubleSpinBoxScale);
        connect(t_pDoubleSpinBoxScale,static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
                this,&QuickControlWidget::updateSpinBoxScaling);
        t_pGridLayout->addWidget(t_pDoubleSpinBoxScale,i+1,0,1,1);

        QSlider* t_pHorizontalSlider = new QSlider(Qt::Horizontal);
        t_pHorizontalSlider->setMinimum(1);
        t_pHorizontalSlider->setMaximum(1000);
        t_pHorizontalSlider->setSingleStep(1);
        t_pHorizontalSlider->setPageStep(1);
        t_pHorizontalSlider->setValue(m_qMapChScaling.value(FIFFV_STIM_CH)/10);
        m_qMapScalingSlider.insert(FIFFV_STIM_CH,t_pHorizontalSlider);
        connect(t_pHorizontalSlider,static_cast<void (QSlider::*)(int)>(&QSlider::valueChanged),
                this,&QuickControlWidget::updateSliderScaling);
        t_pGridLayout->addWidget(t_pHorizontalSlider,i+1,1,1,1);


        i+=2;
    }

    //MISC
    if(m_qMapChScaling.contains(FIFFV_MISC_CH))
    {
        QLabel* t_pLabelModality = new QLabel;
        t_pLabelModality->setText("MISC");
        t_pGridLayout->addWidget(t_pLabelModality,i,0,1,1);

        QDoubleSpinBox* t_pDoubleSpinBoxScale = new QDoubleSpinBox;
        t_pDoubleSpinBoxScale->setMinimum(0.1);
        t_pDoubleSpinBoxScale->setMaximum(10000);
        t_pDoubleSpinBoxScale->setMaximumWidth(100);
        t_pDoubleSpinBoxScale->setSingleStep(0.1);
        t_pDoubleSpinBoxScale->setDecimals(1);
        t_pDoubleSpinBoxScale->setPrefix("+/- ");
        t_pDoubleSpinBoxScale->setValue(m_qMapChScaling.value(FIFFV_MISC_CH));
        m_qMapScalingDoubleSpinBox.insert(FIFFV_MISC_CH,t_pDoubleSpinBoxScale);
        connect(t_pDoubleSpinBoxScale,static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
                this,&QuickControlWidget::updateSpinBoxScaling);
        t_pGridLayout->addWidget(t_pDoubleSpinBoxScale,i+1,0,1,1);

        QSlider* t_pHorizontalSlider = new QSlider(Qt::Horizontal);
        t_pHorizontalSlider->setMinimum(1);
        t_pHorizontalSlider->setMaximum(10000);
        t_pHorizontalSlider->setSingleStep(1);
        t_pHorizontalSlider->setPageStep(1);
        t_pHorizontalSlider->setValue(m_qMapChScaling.value(FIFFV_MISC_CH)/10);
        m_qMapScalingSlider.insert(FIFFV_MISC_CH,t_pHorizontalSlider);
        connect(t_pHorizontalSlider,static_cast<void (QSlider::*)(int)>(&QSlider::valueChanged),
                this,&QuickControlWidget::updateSliderScaling);
        t_pGridLayout->addWidget(t_pHorizontalSlider,i+1,1,1,1);

        i+=2;
    }

    ui->m_groupBox_scaling->setLayout(t_pGridLayout);
}


//*************************************************************************************************************

void QuickControlWidget::createProjectorGroup()
{
    if(m_pFiffInfo)
    {
        m_qListProjCheckBox.clear();
        // Projection Selection
        QGridLayout *topLayout = new QGridLayout;

        bool bAllActivated = true;

        qint32 i=0;

        for(i; i < m_pFiffInfo->projs.size(); ++i)
        {
            QCheckBox* checkBox = new QCheckBox(m_pFiffInfo->projs[i].desc);
            checkBox->setChecked(m_pFiffInfo->projs[i].active);

            if(m_pFiffInfo->projs[i].active == false)
                bAllActivated = false;

            m_qListProjCheckBox.append(checkBox);

            connect(checkBox, static_cast<void (QCheckBox::*)(bool)>(&QCheckBox::clicked),
                    this, &QuickControlWidget::checkProjStatusChanged);

            topLayout->addWidget(checkBox, i, 0); //+2 because we already added two widgets before the first projector check box

//            if(i>m_pFiffInfo->projs.size()/2)
//                topLayout->addWidget(checkBox, i-rowCount, 1); //+2 because we already added two widgets before the first projector check box
//            else {
//                topLayout->addWidget(checkBox, i, 0); //+2 because we already added two widgets before the first projector check box
//                rowCount++;
//            }
        }

        QFrame* line = new QFrame();
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        topLayout->addWidget(line, i+1, 0);

        m_enableDisableProjectors = new QCheckBox("Enable all");
        m_enableDisableProjectors->setChecked(bAllActivated);
        topLayout->addWidget(m_enableDisableProjectors, i+2, 0);
        connect(m_enableDisableProjectors, static_cast<void (QCheckBox::*)(bool)>(&QCheckBox::clicked),
            this, &QuickControlWidget::enableDisableAllProj);

        ui->m_groupBox_projections->setLayout(topLayout);

        //Set default activation to true
        enableDisableAllProj(true);
    }
}


//*************************************************************************************************************

void QuickControlWidget::createOtherGroup()
{
    //Number of visible channels
    connect(ui->m_doubleSpinBox_numberVisibleChannels, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
            this, &QuickControlWidget::zoomChanged);

    //Window size
    connect(ui->m_spinBox_windowSize, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this, &QuickControlWidget::timeWindowChanged);

    //Trigger detection
    connect(ui->m_checkBox_activateTriggerDetection, static_cast<void (QCheckBox::*)(int)>(&QCheckBox::stateChanged),
            this, &QuickControlWidget::realTimeTriggerActiveChanged);

    QMapIterator<QString, QColor> i(m_qMapTriggerColor);
    while(i.hasNext()) {
        i.next();
        ui->m_comboBox_triggerChannels->addItem(i.key());
    }
    connect(ui->m_comboBox_triggerChannels, static_cast<void (QComboBox::*)(const QString&)>(&QComboBox::currentTextChanged),
            this, &QuickControlWidget::realTimeTriggerCurrentChChanged);

    connect(ui->m_pushButton_triggerColor, static_cast<void (QPushButton::*)(bool)>(&QPushButton::clicked),
            this, &QuickControlWidget::realTimeTriggerColorChanged);

    ui->m_pushButton_triggerColor->setAutoFillBackground(true);
    ui->m_pushButton_triggerColor->setFlat(true);
    QPalette* palette1 = new QPalette();
    palette1->setColor(QPalette::Button,QColor(177,0,0));
    ui->m_pushButton_triggerColor->setPalette(*palette1);
    ui->m_pushButton_triggerColor->update();

    connect(ui->m_doubleSpinBox_detectionThresholdFirst, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged),
            this, &QuickControlWidget::realTimeTriggerThresholdChanged);

    connect(ui->m_spinBox_detectionThresholdSecond, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
                this, &QuickControlWidget::realTimeTriggerThresholdChanged);

    connect(ui->m_pushButton_resetNumberTriggers, static_cast<void (QPushButton::*)(bool)>(&QPushButton::clicked),
            this, &QuickControlWidget::onResetTriggerNumbers);

    if(!m_bTriggerDetection)
        ui->m_tabWidget_viewOptions->removeTab(1);

    //opacity
    connect(ui->m_horizontalSlider_opacity, &QSlider::valueChanged,
            this, &QuickControlWidget::onOpacityChange);

    //Distance for timer spacer
    connect(ui->m_comboBox_distaceTimeSpacer, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &QuickControlWidget::onDistanceTimeSpacerChanged);
}


//*************************************************************************************************************

void QuickControlWidget::createModalityGroup()
{
    m_pModalityGroupBox = new QGroupBox;
    m_pModalityGroupBox->setTitle("Modalities");

    m_qListModalities.clear();
    bool hasMag = false;
    bool hasGrad = false;
    bool hasEEG = false;
    bool hasEOG = false;
    bool hasMISC = false;
    for(qint32 i = 0; i < m_pFiffInfo->nchan; ++i)
    {
        if(m_pFiffInfo->chs[i].kind == FIFFV_MEG_CH)
        {
            if(!hasMag && m_pFiffInfo->chs[i].unit == FIFF_UNIT_T)
                hasMag = true;
            else if(!hasGrad &&  m_pFiffInfo->chs[i].unit == FIFF_UNIT_T_M)
                hasGrad = true;
        }
        else if(!hasEEG && m_pFiffInfo->chs[i].kind == FIFFV_EEG_CH)
            hasEEG = true;
        else if(!hasEOG && m_pFiffInfo->chs[i].kind == FIFFV_EOG_CH)
            hasEOG = true;
        else if(!hasMISC && m_pFiffInfo->chs[i].kind == FIFFV_MISC_CH)
            hasMISC = true;
    }

    bool sel = true;
    float val = 1e-11f;

    if(hasMag)
        m_qListModalities.append(Modality("MAG",sel,val));
    if(hasGrad)
        m_qListModalities.append(Modality("GRAD",sel,val));
    if(hasEEG)
        m_qListModalities.append(Modality("EEG",false,val));
    if(hasEOG)
        m_qListModalities.append(Modality("EOG",false,val));
    if(hasMISC)
        m_qListModalities.append(Modality("MISC",false,val));

    QGridLayout* t_pGridLayout = new QGridLayout;

    for(qint32 i = 0; i < m_qListModalities.size(); ++i)
    {
        QString mod = m_qListModalities[i].m_sName;

        QLabel* t_pLabelModality = new QLabel;
        t_pLabelModality->setText(mod);
        t_pGridLayout->addWidget(t_pLabelModality,i,0,1,1);

        QCheckBox* t_pCheckBoxModality = new QCheckBox;
        t_pCheckBoxModality->setChecked(m_qListModalities[i].m_bActive);
        m_qListModalityCheckBox << t_pCheckBoxModality;
        connect(t_pCheckBoxModality,&QCheckBox::stateChanged,
                this,&QuickControlWidget::updateModalityCheckbox);
        t_pGridLayout->addWidget(t_pCheckBoxModality,i,1,1,1);

    }

    m_pModalityGroupBox->setLayout(t_pGridLayout);

    //Decide where to put the group box depending on already available group boxes
    if(!m_bView && m_bFilter)
        ui->m_gridLayout_masterLayout->addWidget(m_pModalityGroupBox, ui->m_gridLayout_masterLayout->rowCount()-1, 0, 1, 1);

    if((m_bView && m_bFilter) || (!m_bView && !m_bFilter))
        ui->m_gridLayout_masterLayout->addWidget(m_pModalityGroupBox, ui->m_gridLayout_masterLayout->rowCount(), 0, 1, 2);

    if(m_bView && !m_bFilter)
        ui->m_gridLayout_masterLayout->addWidget(m_pModalityGroupBox, ui->m_gridLayout_masterLayout->rowCount()-1, 1, 1, 1);
}


//*************************************************************************************************************

void QuickControlWidget::createCompensatorGroup()
{
    if(m_pFiffInfo)
    {
        m_pCompSignalMapper = new QSignalMapper(this);

        m_qListCompCheckBox.clear();

        // Compensation Selection
        QGridLayout *topLayout = new QGridLayout;

        qint32 i=0;

        for(i; i < m_pFiffInfo->comps.size(); ++i)
        {
            QString numStr;
            QCheckBox* checkBox = new QCheckBox(numStr.setNum(m_pFiffInfo->comps[i].kind));

            m_qListCompCheckBox.append(checkBox);

            connect(checkBox, SIGNAL(clicked()),
                        m_pCompSignalMapper, SLOT(map()));

            m_pCompSignalMapper->setMapping(checkBox, numStr);

            topLayout->addWidget(checkBox, i, 0);

        }

        connect(m_pCompSignalMapper, SIGNAL(mapped(const QString &)),
                    this, SIGNAL(compClicked(const QString &)));

        connect(this, &QuickControlWidget::compClicked,
                this, &QuickControlWidget::checkCompStatusChanged);

        ui->m_groupBox_compensators->setLayout(topLayout);
    }
}


//*************************************************************************************************************

void QuickControlWidget::onTimeWindowChanged(int value)
{
    emit timeWindowChanged(value);
}


//*************************************************************************************************************

void QuickControlWidget::onZoomChanged(double value)
{
    emit zoomChanged(value);
}


//*************************************************************************************************************

void QuickControlWidget::checkProjStatusChanged(bool status)
{
    Q_UNUSED(status)

    bool bAllActivated = true;

    for(qint32 i = 0; i < m_qListProjCheckBox.size(); ++i) {
        if(m_qListProjCheckBox[i]->isChecked() == false)
            bAllActivated = false;

        this->m_pFiffInfo->projs[i].active = m_qListProjCheckBox[i]->isChecked();
    }

    m_enableDisableProjectors->setChecked(bAllActivated);

    emit projSelectionChanged();

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::enableDisableAllProj(bool status)
{
    //Set all checkboxes to status
    for(int i=0; i<m_qListProjCheckBox.size(); i++)
        m_qListProjCheckBox.at(i)->setChecked(status);

    //Set all projection activation states to status
    for(int i=0; i < m_pFiffInfo->projs.size(); ++i)
        m_pFiffInfo->projs[i].active = status;

    m_enableDisableProjectors->setChecked(status);

    emit projSelectionChanged();

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::checkCompStatusChanged(const QString & compName)
{
    qDebug()<<compName;

    bool currentState;

    for(int i = 0; i < m_qListCompCheckBox.size(); ++i)
        if(m_qListCompCheckBox[i]->text() != compName)
            m_qListCompCheckBox[i]->setChecked(false);
        else
            currentState = m_qListCompCheckBox[i]->isChecked();

    if(currentState)
        emit compSelectionChanged(compName.toInt());
    else //If none selected
        emit compSelectionChanged(0);

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::updateSpinBoxScaling(double value)
{
    Q_UNUSED(value)

    QMap<qint32, QDoubleSpinBox*>::iterator it;
    for (it = m_qMapScalingDoubleSpinBox.begin(); it != m_qMapScalingDoubleSpinBox.end(); ++it)
    {
        double scaleValue = 0;

        switch(it.key())
        {
            case FIFF_UNIT_T:
                //MAG
                scaleValue = 1e-12;
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*10);
                break;
            case FIFF_UNIT_T_M:
                //GRAD
                scaleValue = 1e-15 * 100; //*100 because data in fiff files is stored as fT/m not fT/cm
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*1);
                break;
            case FIFFV_EEG_CH:
                //EEG
                scaleValue = 1e-06;
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*10);
                break;
            case FIFFV_EOG_CH:
                //EOG
                scaleValue = 1e-06;
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*10);
                break;
            case FIFFV_EMG_CH:
                //EMG
                scaleValue = 1e-03;
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*10);
                break;
            case FIFFV_ECG_CH:
                //ECG
                scaleValue = 1e-03;
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*10);
                break;
            case FIFFV_MISC_CH:
                //MISC
                scaleValue = 1;
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*10);
                break;
            case FIFFV_STIM_CH:
                //STIM
                scaleValue = 1;
                m_qMapScalingSlider[it.key()]->setValue(it.value()->value()*10);
                break;
            default:
                scaleValue = 1.0;
        }

        //if(m_qMapScalingSlider[it.key()]->maximum()<it.value()->value()*10)
            m_qMapChScaling.insert(it.key(), it.value()->value() * scaleValue);
//        qDebug()<<"m_pRTMSAW->m_qMapChScaling[it.key()]" << m_pRTMSAW->m_qMapChScaling[it.key()];
    }

    emit scalingChanged(m_qMapChScaling);

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::updateSliderScaling(int value)
{
    Q_UNUSED(value)

    QMap<qint32, QDoubleSpinBox*>::iterator it;
    for (it = m_qMapScalingDoubleSpinBox.begin(); it != m_qMapScalingDoubleSpinBox.end(); ++it)
    {
        double scaleValue = 0;

        switch(it.key())
        {
            case FIFF_UNIT_T:
                //MAG
                scaleValue = 1e-12;
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/10);
                break;
            case FIFF_UNIT_T_M:
                //GRAD
                scaleValue = 1e-15 * 100; //*100 because data in fiff files is stored as fT/m not fT/cm
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/1);
                break;
            case FIFFV_EEG_CH:
                //EEG
                scaleValue = 1e-06;
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/10);
                break;
            case FIFFV_EOG_CH:
                //EOG
                scaleValue = 1e-06;
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/10);
                break;
            case FIFFV_EMG_CH:
                //EMG
                scaleValue = 1e-03;
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/10);
                break;
            case FIFFV_ECG_CH:
                //ECG
                scaleValue = 1e-03;
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/10);
                break;
            case FIFFV_MISC_CH:
                //MISC
                scaleValue = 1;
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/10);
                break;
            case FIFFV_STIM_CH:
                //STIM
                scaleValue = 1;
                it.value()->setValue((double)m_qMapScalingSlider[it.key()]->value()/10);
                break;
            default:
                scaleValue = 1.0;
        }

//        qDebug()<<"m_pRTMSAW->m_qMapChScaling[it.key()]" << m_pRTMSAW->m_qMapChScaling[it.key()];
    }

//    emit scalingChanged();

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::realTimeTriggerActiveChanged(int state)
{
    Q_UNUSED(state);

    emit triggerInfoChanged(m_qMapTriggerColor, ui->m_checkBox_activateTriggerDetection->isChecked(), ui->m_comboBox_triggerChannels->currentText(), ui->m_doubleSpinBox_detectionThresholdFirst->value()*pow(10, ui->m_spinBox_detectionThresholdSecond->value()));
}


//*************************************************************************************************************

void QuickControlWidget::realTimeTriggerColorChanged(bool state)
{
    Q_UNUSED(state);

    QColor color = QColorDialog::getColor(m_qMapTriggerColor[ui->m_comboBox_triggerChannels->currentText()], this, "Set trigger color");

    //Change color of pushbutton
    QPalette* palette1 = new QPalette();
    palette1->setColor(QPalette::Button,color);
    ui->m_pushButton_triggerColor->setPalette(*palette1);
    ui->m_pushButton_triggerColor->update();

    m_qMapTriggerColor[ui->m_comboBox_triggerChannels->currentText()] = color;

    emit triggerInfoChanged(m_qMapTriggerColor, ui->m_checkBox_activateTriggerDetection->isChecked(), ui->m_comboBox_triggerChannels->currentText(), ui->m_doubleSpinBox_detectionThresholdFirst->value()*pow(10, ui->m_spinBox_detectionThresholdSecond->value()));
}


//*************************************************************************************************************

void QuickControlWidget::realTimeTriggerThresholdChanged(double value)
{
    Q_UNUSED(value);

    emit triggerInfoChanged(m_qMapTriggerColor, ui->m_checkBox_activateTriggerDetection->isChecked(), ui->m_comboBox_triggerChannels->currentText(), ui->m_doubleSpinBox_detectionThresholdFirst->value()*pow(10, ui->m_spinBox_detectionThresholdSecond->value()));
}


//*************************************************************************************************************

void QuickControlWidget::realTimeTriggerCurrentChChanged(const QString &value)
{
    //Change color of pushbutton
    QPalette* palette1 = new QPalette();
    palette1->setColor(QPalette::Button,m_qMapTriggerColor[value]);
    ui->m_pushButton_triggerColor->setPalette(*palette1);
    ui->m_pushButton_triggerColor->update();

    emit triggerInfoChanged(m_qMapTriggerColor, ui->m_checkBox_activateTriggerDetection->isChecked(), ui->m_comboBox_triggerChannels->currentText(), ui->m_doubleSpinBox_detectionThresholdFirst->value()*pow(10, ui->m_spinBox_detectionThresholdSecond->value()));
}


//*************************************************************************************************************

void QuickControlWidget::toggleHideAll(bool state)
{
    if(!state) {
        ui->m_groupBox_projections->hide();
        ui->m_groupBox_filter->hide();
        ui->m_groupBox_scaling->hide();
        ui->m_groupBox_view->hide();
        ui->m_groupBox_compensators->hide();

        if(m_bModalitiy)
            m_pModalityGroupBox->hide();

        ui->m_pushButton_hideAll->setText(QString("Maximize - Quick Control - %1").arg(m_sName));
    }
    else {
        if(m_bProjections)
            ui->m_groupBox_projections->show();

        if(m_bFilter)
            ui->m_groupBox_filter->show();

        if(m_bScaling)
            ui->m_groupBox_scaling->show();

        if(m_bView)
            ui->m_groupBox_view->show();

        if(m_bModalitiy)
            m_pModalityGroupBox->show();

        if(m_bCompensator)
            ui->m_groupBox_compensators->show();

        ui->m_pushButton_hideAll->setText(QString("Minimize - Quick Control - %1").arg(m_sName));
    }

    this->adjustSize();
    this->resize(width(), ui->m_pushButton_hideAll->height()-50);
}


//*************************************************************************************************************

void QuickControlWidget::onShowFilterOptions(bool state)
{
//    if(state)
//        m_pShowFilterOptions->setText("Close filter options");
//    else
//        m_pShowFilterOptions->setText("Open filter options");

//    m_pShowFilterOptions->setChecked(state);

//    emit showFilterOptions(state);

    Q_UNUSED(state);
    emit showFilterOptions(true);
}


//*************************************************************************************************************

void QuickControlWidget::updateModalityCheckbox(qint32 state)
{
    Q_UNUSED(state)

    for(qint32 i = 0; i < m_qListModalityCheckBox.size(); ++i)
    {
        if(m_qListModalityCheckBox[i]->isChecked())
            m_qListModalities[i].m_bActive = true;
        else
            m_qListModalities[i].m_bActive = false;
    }

    emit settingsChanged(m_qListModalities);

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::onOpacityChange(qint32 value)
{
    this->setWindowOpacity(1/(100.0/value));
}


//*************************************************************************************************************

void QuickControlWidget::onDistanceTimeSpacerChanged(qint32 value)
{
    switch(value) {
        case 0:
            emit distanceTimeSpacerChanged(100);
        break;

        case 1:
            emit distanceTimeSpacerChanged(200);
        break;

        case 2:
            emit distanceTimeSpacerChanged(500);
        break;

        case 3:
            emit distanceTimeSpacerChanged(1000);
        break;
    }

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::onResetTriggerNumbers()
{
    ui->m_label_numberDetectedTriggers->setText(QString("0"));

    emit resetTriggerCounter();

    emit updateConnectedView();
}


//*************************************************************************************************************

void QuickControlWidget::userFilterToggled(bool state)
{
    Q_UNUSED(state);
    //qDebug()<<"userFilterToggled";
    emit updateConnectedView();
}






