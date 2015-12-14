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



Eigen::MatrixXd generateProjectionVector( const Eigen::MatrixXd& leftSvdMatrix, const Eigen::RowVectorXi& picks, FIFFLIB::fiff_short_t currentProjectionIndex )
{
    std::cout << currentProjectionIndex << ". vector is being generated" << std::endl;
     const FIFFLIB::fiff_int_t picksSize = picks.cols();
     if( picksSize != leftSvdMatrix.rows() )
         std::cout << "WARNING: different dimensionality of channels and Singular value basis" << std::endl;

    //  TODO: This should work as a RowVector as well... test out later
    Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(1, picksSize);

    for( FIFFLIB::fiff_int_t picksIndex = 0; picksIndex != picksSize; ++picksIndex )
    {
        ret(0, picksIndex) = leftSvdMatrix(picksIndex, currentProjectionIndex );
    }
    std::cout << "Projvector " << currentProjectionIndex << " written" << std::endl;
    return ret;
}



FIFFLIB::FiffNamedMatrix::SDPtr generateFiffProjData( const Eigen::MatrixXd& row, const QStringList& channelNames )
{
    FIFFLIB::FiffNamedMatrix::SDPtr ret = FIFFLIB::FiffNamedMatrix::SDPtr(new FIFFLIB::FiffNamedMatrix());

    ret->nrow = row.rows();
    ret->ncol = row.cols();
    ret->col_names = channelNames;       // TODO: set channel names as in generateFiffProj?
    ret->data = row;

    return ret;

}



//  FUNCTION: GenerateFiffProj
//  Set up a list of vectors to be filtered out
//  Fill up projectionVectors with the wanted vectors.
//  The wanted vectors are the ones of maximum impact in the test run.
//  Maximum impact means they belong to the biggest singular values.
//  Singular Values are ordered by decreasing size.
//  Thus the first vectors of leftSvdMatrix are the wanted vectors.
//
//  TODO:   Use global constant, for #projectionvectors
QList<FIFFLIB::FiffProj> generateFiffProj(const FIFFLIB::FiffRawData& raw, const Eigen::RowVectorXi& picks, const Eigen::MatrixXd& data, FIFFLIB::fiff_short_t numberOfProjections)
{
    if( data.rows() != picks.cols() )
        std::cout << "WARNING: unequal dimensionality of Data matrix and number of picked channels!" << std::endl;
    if( numberOfProjections < 0 )
    {
        std::cout << "WARNING: negative number of projections given" << std::endl;
        numberOfProjections = 0;
    }

    QList<FIFFLIB::FiffProj> ret;   //  initialisation of return variable

    //  Apply Singular Value Decomposition on data, to get the left orthogonal matrix
    //  TODO: check for #projection vectors > #cols(U)
    std::cout << "Computing singular value decomposition..." <<std::flush;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(data, Eigen::ComputeThinU);
    Eigen::MatrixXd leftSvdMatrix = svd.matrixU();
    std::cout << "DONE" << std::endl;
    std::cout << "Sngular value basis of dimensionality " << leftSvdMatrix.rows() << "x" << leftSvdMatrix.cols() << std::endl;


    std::cout << "setting up channel names..." << std::flush;
    const FIFFLIB::fiff_int_t picksSize = picks.cols();
    //  all projection vectors have equal channel names
    //  channelNames will be accessed by generateFiffProjData
    std::cout << "Fetching channel names..." << std::flush;
    QStringList channelNames;
    for( FIFFLIB::fiff_int_t i = 0; i != picksSize; ++i )
    {
        channelNames.push_back( raw.info.ch_names[picks[i]] );
    }
    std::cout << "DONE" << std::endl;


    //Add the #NumberOfProjections first colums of leftSvdMatrix to ret
    std::cout << "Generating the SSP projection file..." << std::flush;
    for( FIFFLIB::fiff_short_t  currentProjectionIndex = 0; currentProjectionIndex != numberOfProjections; ++currentProjectionIndex )
    {
        FIFFLIB::FiffProj retElement;

        retElement.kind = 1;    //  kind of projection-vector = 0 means it is not specified
        retElement.active = false;  //  active-flag can/should be activated manually later on

        // TODO - add impact factor
        retElement.desc = QString("%1. vector - impact = ").arg(currentProjectionIndex + 1);

        Eigen::MatrixXd currentProjectionVector = generateProjectionVector( leftSvdMatrix, picks, currentProjectionIndex);

        retElement.data = generateFiffProjData(currentProjectionVector, channelNames);
        ret.push_back(retElement);
    }
    return ret;
}




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

    /**********************************************************************************************************
     * Read the empty room file, generate the projections
     **********************************************************************************************************/

    Eigen::RowVectorXi picks;
    QList<FIFFLIB::FiffProj> projectionVectors;
    // Delete unnecessary data after this
    {
        //  location of the file that will be read
        QFile t_fileEmptyRoom("./MNE-sample-data/MEG/test_data_ssp/151015_131148_4884471_Empty_room_raw.fif");

        bool in_samples = false;

        //
        //   Setup for reading the raw data
        //
        FIFFLIB::FiffRawData raw(t_fileEmptyRoom);

        float from = 60.0f;
        float to = 120.0f;

        //
        //   Set up pick list: MEG  - bad channels
        //

        //Eigen::RowVectorXi picks = raw.info.pick_types(want_meg, want_eeg, want_stim, include, raw.info.bads);
        //ToDo : Check for bad channels
        std::cout << "Set up pick list..." << std::flush ;
        picks.resize(1, raw.info.chs.size());
        int k = 0;
        for( int i = 0; i != raw.info.chs.size(); ++i )
            if( raw.info.chs[i].unit == 112)
            {
                picks[k] = i;
                ++k;
            }
        picks.conservativeResize(1, k);
        std::cout << "DONE - dimension = " << picks.cols() << std::endl;

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

        std::cout << "Data was read. A matrix of size " << data.rows() << "x" << data.cols() << std::endl;


        //Just for memorizing the code QString test = QString("v%1_%2").arg(0).arg(1);

        // Calculate the projection vectors via SVD and save them in the appropriate structure
        const FIFFLIB::fiff_int_t numberOfProjections = 8;

        projectionVectors = generateFiffProj(raw, picks, data, numberOfProjections);


        t_fileEmptyRoom.close();
    }

    /*************************************************************************************************************************************************
     *  Writing the projection vectors into a file with spontaneous data
     ************************************************************************************************************************************************/

    //  locations of the file to read and where to write the new file
    QFile t_fileIn("./MNE-sample-data/MEG/test_data_ssp/151023_120648_4884471_Spontaneous1_raw.fif");
    QFile t_fileOut("./MNE-sample-data/MEG/test_data_ssp/ssp_output_spontaneous_with_new_projs_raw.fif");

    //
    //   Setup for reading the raw data
    //

    std::cout <<"FiffRawData \n" << std::endl;
    FIFFLIB::FiffRawData rawIn(t_fileIn);
    std::cout <<"FiffRawData finished \n" << std::endl;

    //use old pick list!


    rawIn.proj.resize(0,0);

    rawIn.info.projs = projectionVectors;
    Eigen::RowVectorXd cals;

    Eigen::MatrixXd data;
    Eigen::MatrixXd times;

    std::cout <<"start_writing_raw \n" << std::endl;
    FIFFLIB::FiffStream::SPtr outfid = FIFFLIB::Fiff::start_writing_raw(t_fileOut,rawIn.info, cals/*, picks*/);
    std::cout <<"start_writing_raw finished \n" << std::endl;

    //
    //   Set up the reading parameters
    //
    FIFFLIB::fiff_int_t from = rawIn.first_samp;
    FIFFLIB::fiff_int_t to = rawIn.last_samp;
    //float from = 60.0f;
    //float to = 120.0f;
    float quantum_sec = 10.0f;//read and write in 10 sec junks
    FIFFLIB::fiff_int_t quantum = ceil(quantum_sec*rawIn.info.sfreq);

    //   Read and write all the data
    //
    bool first_buffer = true;

    FIFFLIB::fiff_int_t first, last;

    for(first = from; first < to; first+=quantum)
    {
        last = first+quantum -1;
        if (last > to)
        {
            last = to;
        }

        if (!rawIn.read_raw_segment(data,times,first,last/*,picks*/))
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
















