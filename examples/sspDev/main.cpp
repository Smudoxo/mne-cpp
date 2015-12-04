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
// New Functions
//=============================================================================================================

//=============================================================================================================



Eigen::MatrixXd generateColumn( const Eigen::MatrixXd& leftSvdMatrix, const Eigen::RowVectorXi& picks, FIFFLIB::fiff_int_t vectorSize, FIFFLIB::fiff_short_t currentProjectionIndex )
{
    //  TODO: Sollte auch als reiner Spaltenvektor gehen... klappt dann die conversation später?
    Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(vectorSize, 1);

    FIFFLIB::fiff_int_t picksSize = picks.size();
    if( picksSize != leftSvdMatrix.rows() )
        std::cout << "ERROR: different dimensionality of picks and SVD-basis" << std::endl;     //  TODO: std::except benutzen
    else
    {
        for( FIFFLIB::fiff_int_t picksIndex = 0; picksIndex != picksSize; ++picksIndex )
        {

            ret(picks[picksIndex], 0) = leftSvdMatrix(picksIndex, currentProjectionIndex );
        }

    //    std::cout << std::endl << ret << std::endl;

    //    std::cin.ignore(std::cin.rdbuf()->in_avail());
    //    std::cin.get();
    }

    return ret;
}



FIFFLIB::FiffNamedMatrix::SDPtr generateFiffProjData( const Eigen::MatrixXd& column, const QStringList& channelNames )
{
    FIFFLIB::FiffNamedMatrix::SDPtr ret = FIFFLIB::FiffNamedMatrix::SDPtr(new FIFFLIB::FiffNamedMatrix());

    ret->nrow = column.rows();
    ret->ncol = 1;
    ret->row_names = channelNames;
    ret->col_names.push_back("");         // TODO: set channel names as in generateFiffProj?
    ret->data = column;

    return ret;

}


/*  Old and Useless
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
*/


//  FUNCTION: GenerateFiffProj
//  Set up a list of vectors to be filtered out
//  Fill up projectionVectors with the wanted vectors.
//  The wanted vectors are the ones of maximum impact in the test run.
//  Maximum impact means they belong to the biggest singular values.
//  Singular Values are ordered by decreasing size.
//  Thus the first vectors of leftSvdMatrix are the wanted vectors.
//
//  TODO:   Use global constant, for #projectionvectors
QList<FIFFLIB::FiffProj> generateFiffProj(const FIFFLIB::FiffRawData& raw, const Eigen::RowVectorXi& picks, const Eigen::MatrixXd& data, FIFFLIB::fiff_int_t numberOfProjections)
{
    QList<FIFFLIB::FiffProj> ret;   //  initialisation of return variable

    //  Apply Singular Value Decomposition on data, to get the left orthogonal matrix
    //
    //  TODO: check for #projection vectors > #cols(U)
    std::cout << "Computing singular value decomposition..." <<std::endl;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(data, Eigen::ComputeThinU);
    Eigen::MatrixXd leftSvdMatrix = svd.matrixU();


    //  all projection vectors have equal channel names
    //  channelNames will be accessed by generateFiffProjData
    std::cout << "Fetching channel names..." << std::endl;
    QStringList channelNames = raw.info.ch_names;

    //Add the NumberOfProjections first colums of leftSvdMatrix to ret
    std::cout << "Generating the SSP projection file..." << std::endl;
    for( FIFFLIB::fiff_short_t  currentProjectionIndex = 0; currentProjectionIndex != numberOfProjections; ++currentProjectionIndex )
    {
        FIFFLIB::FiffProj retElement;

        retElement.kind = 0;    //  kind of projection-vector = 0 means it is not specified
        retElement.active = false;  //  active-flag can/should be activated manually later on

        // TODO - add impact factor
        retElement.desc = QString("%1. vector - impact = ").arg(currentProjectionIndex + 1);

        Eigen::MatrixXd currentProjectionColumn = generateColumn( leftSvdMatrix, picks, channelNames.size(), currentProjectionIndex);

        retElement.data = generateFiffProjData(currentProjectionColumn, channelNames);
        ret.push_back(retElement);
    }
    return ret;
}
/*  Useless and Old
Eigen::MatrixXd generateMatrixfromFiffProj( const QList<FIFFLIB::FiffProj>& Proj )
{
    Eigen::MatrixXd ret;
    for (int j = 0; j != Proj.size(); ++j )
    {
        Eigen::MatrixXd temp = (Proj.at(j).data)->data;
        for( int i = 0; i != temp.cols(); ++i )
        {
            ret(i, j) = temp(0,i);
        }
    }
    return ret;
}
*/




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

    //
    //  location of file to read
    //
    QFile t_fileRaw("./MNE-sample-data/MEG/sample/sample_audvis_raw.fif");

    bool in_samples = true;

    //
    //   Setup for reading the raw data
    //
    FIFFLIB::FiffRawData raw(t_fileRaw);

    FIFFLIB::fiff_int_t from = raw.first_samp;
    FIFFLIB::fiff_int_t to = raw.last_samp;

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

    //  Print a message if reading wasn't successful
    if (!readSuccessful)
    {
        printf("Could not read raw segment.\n");
        return -1;
    }

    printf("Read %d samples.\n",(qint32)data.cols());


    //Print out a sample of the data,
    std::cout << data.block(0,0,5,5) << std::endl;


    //Just for memorizing the code QString test = QString("v%1_%2").arg(0).arg(1);


    // Calculate the projection vectors via SVD and save them in the appropriate structure
    const FIFFLIB::fiff_int_t numberOfProjections = 3;
    QList<FIFFLIB::FiffProj> projectionVectors = generateFiffProj(raw, picks, data, numberOfProjections);


    //raw.proj = generateMatrixfromFiffProj(projectionVectors);
    //raw.info.projs = projectionVectors;


    /****************************************************************************************************
     ****************************************************************************************************
     * Writing the new, raw file, for testing reasons
     ****************************************************************************************************
     ****************************************************************************************************
     */
        QFile t_fileOut("./MNE-sample-data/MEG/sample/test_output_with_new_projs.fif");


        Eigen::RowVectorXd cals;

        FIFFLIB::FiffStream::SPtr outfid = FIFFLIB::Fiff::start_writing_raw(t_fileOut,raw.info, cals/*, picks*/);


        float quantum_sec = 10.0f;//read and write in 10 sec junks
        FIFFLIB::fiff_int_t quantum = ceil(quantum_sec*raw.info.sfreq);
        //
        //   To read the whole file at once set
        //
        //quantum     = to - from + 1;
        //
        //
        //   Read and write all the data
        //
        bool first_buffer = true;

        FIFFLIB::fiff_int_t first, last;


        for(first = from; first < to; first+=quantum)
        {
            last = first+quantum-1;
            if (last > to)
            {
                last = to;
            }

            if (!raw.read_raw_segment(data,times,first,last/*,picks*/))
            {
                    printf("error during read_raw_segment\n");
                    return -1;
            }
            //
            //   You can add your own miracle here
            //
            printf("Writing...");
            if (first_buffer)
            {
               if (first > 0)
                   outfid->write_int(FIFF_FIRST_SAMPLE,&first);
               first_buffer = false;
            }
            outfid->write_raw_buffer(data,cals);
            printf("[done]\n");
        }

        outfid->finish_writing_raw();

        printf("Finished\n");

        return 0;//a.exec();
}
















