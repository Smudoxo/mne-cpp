//=============================================================================================================
/**
* @file     rawdelegate.cpp
* @author   Florian Schlembach <florian.schlembach@tu-ilmenau.de>;
*           Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>;
*           Jens Haueisen <jens.haueisen@tu-ilmenau.de>
* @version  1.0
* @date     January, 2014
*
* @section  LICENSE
*
* Copyright (C) 2014, Florian Schlembach, Christoph Dinh and Matti Hamalainen. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that
* the following conditions are met:
*     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
*       following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
*       the following disclaimer in the documentation and/or other materials provided with the distribution.
*     * Neither the name of the Massachusetts General Hospital nor the names of its contributors may be used
*       to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
* PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL MASSACHUSETTS GENERAL HOSPITAL BE LIABLE FOR ANY DIRECT,
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*
* @brief    Implementation of delegate of mne_browse_raw_qt
*
*/


//=============================================================================================================
// INCLUDES

#include "rawdelegate.h"
#include "types_settings.h"

//Qt
#include <QPointF>
#include <QRect>

//=============================================================================================================
// NAMESPACES

using namespace MNE_BROWSE_RAW_QT;
using namespace Eigen;
using namespace MNELIB;

//*************************************************************************************************************

RawDelegate::RawDelegate(QObject *parent)
: m_dPlotHeight(70)
, m_dDx(1)
, m_nhlines(6)
{
}

//*************************************************************************************************************

void RawDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    switch(index.column()) {
    case 0: //chnames
        painter->save();

//        qDebug() << "option.rect.x" << option.rect.x() << "y" << option.rect.y() << "w" << option.rect.width() << "h" << option.rect.height();
        painter->rotate(-90);
        painter->drawText(QRectF(-option.rect.y()-m_dPlotHeight,0,m_dPlotHeight,20),Qt::AlignCenter,index.model()->data(index,Qt::DisplayRole).toString());
        painter->restore();
        break;
    case 1: //data plot
        painter->save();

        //Get data
        QVariant variant = index.model()->data(index,Qt::DisplayRole);
        QList<RowVectorPair> listPairs = variant.value<QList<RowVectorPair> >();
        const RawModel* t_rawModel = (static_cast<const RawModel*>(index.model()));

        QPainterPath path(QPointF(option.rect.x()+t_rawModel->relFiffCursor()-1,option.rect.y()));

        //Plot grid
        painter->setRenderHint(QPainter::Antialiasing, false);
        createGridPath(path,listPairs);

        painter->save();
        QPen pen;
        pen.setStyle(Qt::DotLine);
        pen.setWidthF(0.5);
        painter->setPen(pen);
        painter->drawPath(path);
        painter->restore();

        //Plot data path
        path = QPainterPath(QPointF(option.rect.x()+t_rawModel->relFiffCursor(),option.rect.y()));
        createPlotPath(index,path,listPairs);

        painter->translate(0,m_dPlotHeight/2);

        painter->setRenderHint(QPainter::Antialiasing, true);
        painter->drawPath(path);

        painter->restore();

        break;
    }

}

//*************************************************************************************************************

QSize RawDelegate::sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    QSize size;

    switch(index.column()) {
    case 0:
        size = QSize(20,m_dPlotHeight);
        break;
    case 1:
        QList<RowVectorPair> listPairs = index.model()->data(index).value<QList<RowVectorPair> >();
        qint32 nsamples = (static_cast<const RawModel*>(index.model()))->lastSample()-(static_cast<const RawModel*>(index.model()))->firstSample();

        size = QSize(nsamples*m_dDx,m_dPlotHeight);
        break;
    }

    return size;
}

//=============================================================================================================

void RawDelegate::createPlotPath(const QModelIndex &index, QPainterPath& path, QList<RowVectorPair>& listPairs) const
{
    //get maximum range of respective channel type (range value in FiffChInfo does not seem to contain a reasonable value)
    qint32 kind = (static_cast<const RawModel*>(index.model()))->m_chinfolist[index.row()].kind;
    double dMaxValue = 1e-9;

    switch(kind) {
    case FIFFV_MEG_CH: {
        qint32 unit = (static_cast<const RawModel*>(index.model()))->m_pfiffIO->m_qlistRaw[0]->info.chs[index.row()].unit;
        if(unit == FIFF_UNIT_T_M) {
            dMaxValue = MAX_MEG_UNIT_T_M;
        }
        else if(unit == FIFF_UNIT_T)
            dMaxValue = MAX_MEG_UNIT_T;
        break;
    }
    case FIFFV_EEG_CH: {
        dMaxValue = MAX_EEG;
        break;
    }
    case FIFFV_EOG_CH: {
        dMaxValue = MAX_EOG;
        break;
    }
    case FIFFV_STIM_CH: {
        dMaxValue = MAX_STIM;
        break;
    }
    }

    double dValue;
//    double dMaxValue = (static_cast<const RawModel*>(index.model()))->maxDataValue(index.row());
    double dScaleY = m_dPlotHeight/(2*dMaxValue);

    double y_base = path.currentPosition().y();
    QPointF qSamplePosition;

    //plot all rows from list of pairs
    for(qint8 i=0; i < listPairs.size(); ++i) {
        //create lines from one to the next sample
        for(qint32 j=0; j < listPairs[i].second; ++j) //ToDo: check whether -1 is necessary
        {
            double val = *(listPairs[i].first+j);
            dValue = val*dScaleY;

            double newY = y_base+dValue;

            qSamplePosition.setY(newY);
            qSamplePosition.setX(path.currentPosition().x()+m_dDx);

            path.lineTo(qSamplePosition);
        }
    }

//    qDebug("Plot-PainterPath created!");
}

//*************************************************************************************************************

void RawDelegate::createGridPath(QPainterPath& path, QList<RowVectorPair>& listPairs) const
{
    //horizontal lines
    double distance = m_dPlotHeight/m_nhlines;

    QPointF startpos = path.currentPosition();
    QPointF endpoint(path.currentPosition().x()+listPairs[0].second*listPairs.size()*m_dDx,path.currentPosition().y());

    for(qint8 i=0; i < m_nhlines-1; ++i) {
        endpoint.setY(endpoint.y()+distance);
        path.moveTo(startpos.x(),endpoint.y());
        path.lineTo(endpoint);
    }

//    qDebug("Grid-PainterPath created!");
}
