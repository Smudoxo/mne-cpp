//=============================================================================================================
/**
* @file     main.cpp
* @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
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
* @brief    Example of reading raw data
*
*/


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <iostream>
#include <vector>
#include <math.h>

#include <fiff/fiff.h>


//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtCore/QCoreApplication>


//*************************************************************************************************************
//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>
#include <Eigen/SVD>


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

//using namespace FIFFLIB;
//using namespace Eigen;


//*************************************************************************************************************
//=============================================================================================================
// MAIN
//=============================================================================================================

//=============================================================================================================
/**
* The function main marks the entry point of the program.
* By default, main has the storage class extern.
*
* @param [in] argc (argument count) is an integer that indicates how many arguments were entered on the command line when the program was started.
* @param [in] argv (argument vector) is an array of pointers to arrays of character objects. The array objects are null-terminated strings, representing the arguments that were entered on the command line when the program was started.
* @return the value that was set to exit() (which is 0 if exit() is called via quit()).
*/
int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    //Warum MatrixXd - dynamische Größe ist nicht notwendig...

    /*
     * Unklarheit, Ausgabe der letzten Zeile fehlt
     *
    Eigen::MatrixXd testMatrix;
    testMatrix.resize(3,3);
    testMatrix.setZero();

    testMatrix(0,0) = 3;
    testMatrix(1,1) = 1;
    testMatrix(2,2) = 5;

    std::cout << testMatrix;
    */

    //
    //  location of file to read
    //
    QFile t_fileRaw("./MNE-sample-data/MEG/sample/sample_audvis_raw.fif");

    float from = 42.956f;
    float to = 320.670f;

    bool in_samples = false;


    //
    //   Setup for reading the raw data
    //
    FIFFLIB::FiffRawData raw(t_fileRaw);

    //
    //   Set up pick list: MEG  - bad channels
    //

    QStringList include;  //empty additnional includes
    bool want_meg   = true;
    bool want_eeg   = false;
    bool want_stim  = false;

    Eigen::RowVectorXi picks = raw.info.pick_types(want_meg, want_eeg, want_stim, include, raw.info.bads);

    //
    //  delete projections of raw file
    //  otherwise the projection will be carried out while reading ot the file
    //
    raw.proj.resize(0,0);

    //
    //   Read a data segment
    //   times output argument is optional
    //
    bool readSuccessful = false;
    Eigen::MatrixXd data;
    Eigen::MatrixXd times;
    if (in_samples)
        readSuccessful = raw.read_raw_segment(data, times, (qint32)from, (qint32)to, picks);
    else
        readSuccessful = raw.read_raw_segment_times(data, times, from, to, picks);

    if (!readSuccessful)
    {
        printf("Could not read raw segment.\n");
        return -1;
    }

    printf("Read %d samples.\n",(qint32)data.cols());


    //Print out a sample of the data,
    std::cout << data.block(0,0,5,5) << std::endl;


    //DELETE QString test = QString("v%1_%2").arg(0).arg(1);


    //  Apply Singular Value Decomposition on data, to get the left orthogonal matrix
    //
    //  TODO: check for #projection vectors > #cols(U)

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(data, Eigen::ComputeThinU);
    Eigen::MatrixXd leftSvdMatrixU = svd.matrixU();

    qDebug()<<"SVD calculated";
    //  FUNCTION: GenerateFiffProj
    //  Set up a list of vectors to be filtered out
    //  Fill up projectionVectors with the wanted vectors.
    //  The wanted vectors are the ones of maximum impact in the test run.
    //  Maximum impact means they belong to the biggest singular values.
    //  Singular Values are ordered by decreasing size.
    //  Thus the first vectors of leftSvdMatrix are the wanted vectors.
    //
    //  TODO:   Use global constant, for #projectionvectors
    //  Wie konstruiere ich einen FiffProj aus einer MatrixXd?

    QList<FIFFLIB::FiffProj> projectionVectors;
    //const FIFFLIB::fiff_int_t numberOfProjections = 10;




    return a.exec();
}



Eigen::MatrixXd generateColumn( Eigen::MatrixXd leftSvdMatrix, FIFFLIB::fiff_short_t currentProjectionIndex )
{
    const FIFFLIB::fiff_int_t rowSize = leftSvdMatrix.rows();
    Eigen::MatrixXd ret(rowSize, 1);
    for( FIFFLIB::fiff_int_t rowIndex = 0; rowIndex != rowSize; ++rowIndex )
    {
        ret(rowIndex, 0) = leftSvdMatrix(rowIndex, currentProjectionIndex );
    }

    return ret;
}



FIFFLIB::FiffNamedMatrix::SDPtr generateFiffProjData(  const Eigen::RowVectorXi& picks,
                                                const Eigen::MatrixXd& column,
                                                const QStringList& channelNames )
{
    FIFFLIB::FiffNamedMatrix::SDPtr ret = FIFFLIB::FiffNamedMatrix::SDPtr(new FIFFLIB::FiffNamedMatrix());

    ret->nrow = picks.cols();    //  since picks contains all the used indices of the raw file
    ret->ncol = 1;
    ret->row_names = channelNames;
    ret->col_names.push_back("");         // TODO: set channel names as in generateFiffProj?
    ret->data = column;

    return ret;

}



QStringList generateChannelnames(const FIFFLIB::FiffRawData& raw, const Eigen::RowVectorXi& picks)
{
    QStringList ret;
    for( FIFFLIB::fiff_int_t i = 0; i != picks.cols(); ++i )
    {
        if(picks(i) < raw.info.ch_names.size())
            ret.push_back( raw.info.ch_names.at(picks(i)));
    }
    return ret;
}



QList<FIFFLIB::FiffProj> generateFiffProj(const FIFFLIB::FiffRawData& raw, const Eigen::RowVectorXi& picks, const Eigen::MatrixXd& data, FIFFLIB::fiff_int_t numberOfProjections)
{
    QList<FIFFLIB::FiffProj> ret;   //  initialisation of return variable

    //  Apply Singular Value Decomposition on data, to get the left orthogonal matrix
    //
    //  TODO: check for #projection vectors > #cols(U)
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(data, Eigen::ComputeThinU);
    Eigen::MatrixXd leftSvdMatrix = svd.matrixU();

    //  all projection vectors have equal channel names
    //  channelNames will be accessed by generateFiffProjData
    QStringList channelNames = generateChannelnames(raw, picks);

    //Add the NumberOfProjections first colums of leftSvdMatrix to ret
    for( FIFFLIB::fiff_short_t  currentProjectionIndex = 0;
                                currentProjectionIndex != numberOfProjections;
                                ++currentProjectionIndex )
    {
        FIFFLIB::FiffProj retElement;

        retElement.kind = 0;    //  kind of projection-vector = 0 means it is not specified

        retElement.active = false;  //  active-flag can/should be activated manually later on

        // TODO - add impact factor
        retElement.desc = QString("%1. vector - impact = ").arg(currentProjectionIndex + 1);

        Eigen::MatrixXd currentProjectionColumn = generateColumn( leftSvdMatrix, currentProjectionIndex);

        retElement.data = generateFiffProjData(picks, currentProjectionColumn, channelNames);
        ret.push_back(retElement);
    }
    return ret;
}
















