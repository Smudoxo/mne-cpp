//=============================================================================================================
/**
* @file     fiff_info.cpp
* @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
*           Matti H�m�l�inen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     July, 2012
*
* @section  LICENSE
*
* Copyright (C) 2012, Christoph Dinh and Matti Hamalainen. All rights reserved.
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
* @brief    Contains the implementation of the FiffInfo Class.
*
*/

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "fiff_info.h"


//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================

#include "../3rdParty/Eigen/SVD"


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace FIFFLIB;
using namespace Eigen;


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

FiffInfo::FiffInfo()
: acq_pars("")
, acq_stim("")
, filename("")
{
    meas_date[0] = -1;
}


//*************************************************************************************************************

FiffInfo::~FiffInfo()
{
    qint32 i;
    for (i = 0; i < projs.size(); ++i)
        if(projs[i])
            delete projs[i];
    for (i = 0; i < comps.size(); ++i)
        if(comps[i])
            delete comps[i];
}


//*************************************************************************************************************

qint32 FiffInfo::get_current_comp()
{
    qint32 comp = 0;
    qint32 first_comp = -1;

    qint32 k = 0;
    for (k = 0; k < this->nchan; ++k)
    {
        if (this->chs[k].kind == FIFFV_MEG_CH)
        {
            comp = this->chs[k].coil_type >> 16;
            if (first_comp < 0)
                first_comp = comp;
            else if (comp != first_comp)
                printf("Compensation is not set equally on all MEG channels");
        }
    }
    return comp;
}


//*************************************************************************************************************

bool FiffInfo::make_compensator(fiff_int_t from, fiff_int_t to, FiffCtfComp& ctf_comp, bool exclude_comp_chs)
{
    qDebug() << "make_compensator not debugged jet";

    MatrixXf C1, C2, comp_tmp;

    qDebug() << "Todo add all need ctf variables.";
    if(ctf_comp.data)
        delete ctf_comp.data;
    ctf_comp.data = new FiffNamedMatrix();

    if (from == to)
    {
        ctf_comp.data->data = MatrixXf::Zero(this->nchan, this->nchan);
        return false;
    }

    if (from == 0)
        C1 = MatrixXf::Zero(this->nchan,this->nchan);
    else
    {
        if (!this->make_compensator(from, C1))
        {
            printf("Cannot create compensator C1\n");
            printf("Desired compensation matrix (kind = %d) not found\n",from);
            return false;
        }
    }

    if (to == 0)
        C2 = MatrixXf::Zero(this->nchan,this->nchan);
    else
    {
        if (!this->make_compensator(to, C2))
        {
            printf("Cannot create compensator C2\n");
            printf("Desired compensation matrix (kind = %d) not found\n",to);
            return false;
        }
    }
    //
    //   s_orig = s_from + C1*s_from = (I + C1)*s_from
    //   s_to   = s_orig - C2*s_orig = (I - C2)*s_orig
    //   s_to   = (I - C2)*(I + C1)*s_from = (I + C1 - C2 - C2*C1)*s_from
    //
    comp_tmp = MatrixXf::Identity(this->nchan,this->nchan) + C1 - C2 - C2*C1;

    qint32 k;
    if (exclude_comp_chs)
    {
        VectorXi pick  = MatrixXi::Zero(1,this->nchan);
        qint32 npick = 0;
        for (k = 0; k < this->nchan; ++k)
        {
            if (this->chs[k].kind != FIFFV_REF_MEG_CH)
            {
                pick(npick) = k;
                ++npick;
            }
        }
        if (npick == 0)
        {
            printf("Nothing remains after excluding the compensation channels\n");
            return false;
        }

        ctf_comp.data->data.resize(npick,this->nchan);
        for (k = 0; k < npick; ++k)
            ctf_comp.data->data.row(k) = comp_tmp.block(pick(k), 0, 1, this->nchan);
    }
    return true;
}


//*************************************************************************************************************

bool FiffInfo::make_compensator(fiff_int_t kind, MatrixXf& this_comp)//private method
{
    qDebug() << "make_compensator not debugged jet";
    FiffNamedMatrix* this_data;
    MatrixXf presel, postsel;
    qint32 k, col, c, ch, row, row_ch, channelAvailable;
    for (k = 0; k < this->comps.size(); ++k)
    {
        if (this->comps[k]->kind == kind)
        {                this_data = this->comps[k]->data;
            //
            //   Create the preselector
            //
            presel  = MatrixXf::Zero(this_data->ncol,this->nchan);
            for(col = 0; col < this_data->ncol; ++col)
            {
                channelAvailable = 0;
                for (c = 0; c < this->ch_names.size(); ++c)
                {
                    if (QString::compare(this_data->col_names.at(col),this->ch_names.at(c)) == 0)
                    {
                        ++channelAvailable;
                        ch = c;
                    }
                }
                if (channelAvailable == 0)
                {
                    printf("Channel %s is not available in data\n",this_data->col_names.at(col).toUtf8().constData());
                    return false;
                }
                else if (channelAvailable > 1)
                {
                    printf("Ambiguous channel %s",this_data->col_names.at(col).toUtf8().constData());
                    return false;
                }
                presel(col,ch) = 1.0;
            }
            //
            //   Create the postselector
            //
            postsel = MatrixXf::Zero(this->nchan,this_data->nrow);
            for (c = 0; c  < this->nchan; ++c)
            {
                channelAvailable = 0;
                for (row = 0; row < this_data->row_names.size(); ++row)
                {
                    if (QString::compare(this_data->col_names.at(c),this->ch_names.at(row)) == 0)
                    {
                        ++channelAvailable;
                        row_ch = row;
                    }
                }
                if (channelAvailable > 1)
                {
                    printf("Ambiguous channel %s", this->ch_names.at(c).toUtf8().constData());
                    return false;
                }
                else if (channelAvailable == 1)
                {
                    postsel(c,row_ch) = 1.0;
                }
            }
            this_comp = postsel*this_data->data*presel;
            return true;
        }
    }
    this_comp = defaultMatrixXf;
    return false;
}


//*************************************************************************************************************

fiff_int_t FiffInfo::make_projector(QList<FiffProj*>& projs, QStringList& ch_names, MatrixXf& proj, QStringList& bads, MatrixXf& U)
{
    fiff_int_t nchan = ch_names.size();
    if (nchan == 0)
    {
        printf("No channel names specified\n");
        return 0;
    }

    proj = MatrixXf::Identity(nchan,nchan);
    fiff_int_t nproj = 0;

    //
    //   Check trivial cases first
    //
    if (projs.size() == 0)
        return 0;

    fiff_int_t nactive = 0;
    fiff_int_t nvec    = 0;
    fiff_int_t k, l;
    for (k = 0; k < projs.size(); ++k)
    {
        if (projs[k]->active)
        {
            ++nactive;
            nvec += projs[k]->data->nrow;
        }
    }

    if (nactive == 0)
        return 0;

    //
    //   Pick the appropriate entries
    //
    MatrixXf vecs = MatrixXf::Zero(nchan,nvec);
    nvec = 0;
    fiff_int_t nonzero = 0;
    qint32 p, c, i, j, v;
    double onesize;
    bool isBad = false;
    MatrixXi sel(1, nchan);
    MatrixXi vecSel(1, nchan);
    sel.setConstant(-1);
    vecSel.setConstant(-1);
    for (k = 0; k < projs.size(); ++k)
    {
        sel.resize(1, nchan);
        vecSel.resize(1, nchan);
        sel.setConstant(-1);
        vecSel.setConstant(-1);
        if (projs[k]->active)
        {
            FiffProj* one = projs[k];

            QMap<QString, int> uniqueMap;
            for(l = 0; l < one->data->col_names.size(); ++l)
                uniqueMap[one->data->col_names[l] ] = 0;

            if (one->data->col_names.size() != uniqueMap.keys().size())
            {
                printf("Channel name list in projection item %d contains duplicate items");
                return 0;
            }

            //
            // Get the two selection vectors to pick correct elements from
            // the projection vectors omitting bad channels
            //
            p = 0;
            for (c = 0; c < nchan; ++c)
            {
                for (i = 0; i < one->data->col_names.size(); ++i)
                {
                    if (QString::compare(ch_names.at(c),one->data->col_names.at(i)) == 0)
                    {
                        isBad = false;
                        for (j = 0; j < bads.size(); ++j)
                        {
                            if (QString::compare(ch_names.at(c),bads.at(j)) == 0)
                            {
                                isBad = true;
                            }
                        }

                        if (!isBad && sel(0,p) != c)
                        {
                            sel(0,p) = c;
                            vecSel(0, p) = i;
                            ++p;
                        }

                    }
                }
            }
            sel.conservativeResize(1, p);
            vecSel.conservativeResize(1, p);
            //
            // If there is something to pick, pickit
            //
            if (sel.cols() > 0)
            {
                for (v = 0; v < one->data->nrow; ++v)
                    for (i = 0; i < p; ++i)
                        vecs(sel(0,i),nvec+v) = one->data->data(v,vecSel(i));
                //
                //   Rescale for more straightforward detection of small singular values
                //

                for (v = 0; v < one->data->nrow; ++v)
                {
                    onesize = sqrt((vecs.col(nvec+v).transpose()*vecs.col(nvec+v))(0,0));
                    if (onesize > 0.0)
                    {
                        vecs.col(nvec+v) = vecs.col(nvec+v)/onesize;
                        ++nonzero;
                    }
                }
                nvec += one->data->nrow;
            }
        }
    }
    //
    //   Check whether all of the vectors are exactly zero
    //
    if (nonzero == 0)
        return 0;

    //
    //   Reorthogonalize the vectors
    //
    qDebug() << "Attention Jacobi SVD is used, not the MATLAB lapack version. Since the SVD is not unique the results might be a bit different!";

    JacobiSVD<MatrixXf> svd(vecs.block(0,0,vecs.rows(),nvec), ComputeThinU);

    VectorXf S = svd.singularValues();

    //
    //   Throw away the linearly dependent guys
    //
    for(k = 0; k < nvec; ++k)
    {
        if (S(k)/S(0) < 1e-2)
        {
            nvec = k+1;
            break;
        }
    }

    U = svd.matrixU().block(0, 0, vecs.rows(), nvec);

    //
    //   Here is the celebrated result
    //
    proj -= U*U.transpose();
    nproj = nvec;

    return nproj;
}


//*************************************************************************************************************

QList<FiffChInfo> FiffInfo::set_current_comp(QList<FiffChInfo>& chs, fiff_int_t value)
{
    QList<FiffChInfo> new_chs;
    qint32 k;
    fiff_int_t coil_type;

    for(k = 0; k < chs.size(); ++k)
    {
        new_chs.append(FiffChInfo(&chs[k]));
    }

    qint32 lower_half = 65535;// hex2dec('FFFF');
    for (k = 0; k < chs.size(); ++k)
    {
        if (chs[k].kind == FIFFV_MEG_CH)
        {
            coil_type = chs[k].coil_type & lower_half;
            new_chs[k].coil_type = (coil_type | (value << 16));
        }
    }
    return new_chs;
}