//=============================================================================================================
/**
* @file     brainobject.cpp
* @author   Lorenz Esch <Lorenz.Esch@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     November, 2015
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
* @brief    BrainObject class definition.
*
*/

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "brainobject.h"


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace DISP3DNEWLIB;


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

BrainObject::BrainObject(const Surface &tSurface, const Annotation &tAnnotation, Qt3DCore::QEntity *parent)
: Renderable3DEntity(tSurface.rr(), tSurface.nn(), tSurface.tris(), -tSurface.offset(), parent)
, m_sFilePath(tSurface.filePath())
, m_sSurfFileName(tSurface.fileName())
, m_iHemi(tSurface.hemi())
, m_sSurf(tSurface.surf())
, m_vecCurv(tSurface.curv())
, m_vecOffset(tSurface.offset())
, m_ColorGyri(QColor(125,125,125))
, m_ColorSulci(QColor(50,50,50))
, m_matVert(tSurface.rr())
, m_matTris(tSurface.tris())
, m_matNorm(tSurface.nn())
, m_sAnnotFilePath(tAnnotation.fileName())
{
    //Create color from curvature information
    m_matColorsOrig.resize(m_matVert.rows(), m_matVert.cols());

    for(int i = 0; i<m_matVert.rows() ; i++) {
        if(m_vecCurv[i] >= 0) {
            m_matColorsOrig(i, 0) = m_ColorSulci.redF();
            m_matColorsOrig(i, 1) = m_ColorSulci.greenF();
            m_matColorsOrig(i, 2) = m_ColorSulci.blueF();
        } else {
            m_matColorsOrig(i, 0) = m_ColorGyri.redF();
            m_matColorsOrig(i, 1) = m_ColorGyri.greenF();
            m_matColorsOrig(i, 2) = m_ColorGyri.blueF();
        }
    }

    //Create color from annotation data if annotation is not empty
    if(tAnnotation.getVertices().rows() != 0) {
        tAnnotation.toLabels(tSurface, m_qListLabels, m_qListLabelRGBAs);

        m_matColorsAnnot.resize(m_matVert.rows(), m_matVert.cols());

        for(int i = 0; i<m_qListLabels.size(); i++) {
            FSLIB::Label label = m_qListLabels.at(i);
            for(int j = 0; j<label.vertices.rows(); j++) {
                m_matColorsAnnot(label.vertices(j), 0) = m_qListLabelRGBAs.at(i)(0)/255.0;
                m_matColorsAnnot(label.vertices(j), 1) = m_qListLabelRGBAs.at(i)(1)/255.0;
                m_matColorsAnnot(label.vertices(j), 2) = m_qListLabelRGBAs.at(i)(2)/255.0;
            }
        }
    }
}


//*************************************************************************************************************

BrainObject::~BrainObject()
{
}


//*************************************************************************************************************

void BrainObject::showAnnotation(bool flag)
{
    if(flag && m_matColorsAnnot.rows()!=0)
        this->updateVertColors(m_matColorsAnnot);
    else
        this->updateVertColors(m_matColorsOrig);
}



