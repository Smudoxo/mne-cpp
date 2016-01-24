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
#include <algorithm>

#include <string>

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
// 2015 Version of SSP Generation-function
//=============================================================================================================

//=============================================================================================================



Eigen::MatrixXd generateProjectionVector( const Eigen::MatrixXd& leftSvdMatrix, const Eigen::RowVectorXi& picks, FIFFLIB::fiff_short_t currentProjectionIndex )
{
     const FIFFLIB::fiff_int_t picksSize = picks.cols();
     if( picksSize != leftSvdMatrix.rows() )
         std::cout << "WARNING: different dimensionality of channels and Singular value basis" << std::endl;

    //  TODO: This should work as a RowVector as well... test out later
    Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(1, picksSize);

    for( FIFFLIB::fiff_int_t picksIndex = 0; picksIndex != picksSize; ++picksIndex )
    {
        ret(0, picksIndex) = leftSvdMatrix(picksIndex, currentProjectionIndex );
    }
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
    std::cout << "Singular value basis of dimensionality " << leftSvdMatrix.rows() << "x" << leftSvdMatrix.cols() << std::endl;
    std::cout << "Singular Values are: "
              << svd.singularValues() << std::endl;
    int k;
    std::cin >> k ;


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


/**********************************************************************************************
 *  New Version!
 *  January 2016
 *********************************************************************************************/


/*
 *  CREATE CHANNEL SELECTION
 *
 *  Check each channel, if its kind == 'ch_class'  //   ToDo: Check if kind is right here
 *  and if its name is in the bads
 */
//  ToDo Check that chs contains the necessary data,
//  esp. chs != empty

QStringList createChannelSelection_Matti( const FIFFLIB::FiffCov& cov, FIFFLIB::fiff_int_t ch_class)
{
    QStringList names;
    //  write the names of all channels, which are of the right channel_class and not bad channels, int to 'names'
    for (FIFFLIB::fiff_int_t k = 0; k < cov.names.size(); ++k)
    {
        if (cov.ch_class[k] == ch_class)
        {
            //  ToDo: for speed increase - sort bads and use binary search
            if (!cov.bads.contains(cov.names[k]))
            {
                names.append(cov.names[k]);
            }
        }
    }
    return names;
}

//  ToDo: Ask Lorenz, if this is necessary/desired, to overload this function once more (if covariance matrix is created anyway, this is not necessary
//  Last update at 7. of January 2016
Eigen::RowVectorXi createChannelSelection_Thomas( const QList<FIFFLIB::FiffChInfo>& chs, const QStringList& bads, const FIFFLIB::fiff_int_t& ch_class)
{
    Eigen::RowVectorXi picks = Eigen::RowVectorXi(chs.size());
    int k = 0;
    //  ToDo: Speed increase - sort bad, use binary_search to check for bads
    //  ToDo: switch check for chs[i].unit with check for chs[i].kind, if this is what kind describes
    for( int i = 0; i != chs.size(); ++i )
        //  ToDo: Check if condition <=> chs[i].unit == 112 && ...
        if( chs[i].kind == ch_class && !bads.contains(chs[i].ch_name) )
        {
            picks[k] = i;
            ++k;
        }
    picks.conservativeResize(1, k);

    return picks;
}


Eigen::RowVectorXi createChannelSelection_Thomas( const FIFFLIB::FiffCov& cov, FIFFLIB::fiff_int_t ch_class)
{
    Eigen::RowVectorXi picks = Eigen::RowVectorXi(cov.names.size());
    int k = 0;
    //  ToDo: Speed increase - sort bad, use binary_search to check for bads
    //  ToDo: switch check for chs[i].unit with check for chs[i].kind, if this is what kind describes
    for( int i = 0; i != cov.names.size(); ++i )
        //  ToDo: Check if condition <=> chs[i].unit == 112 && ...
        if( cov.ch_class[i] == ch_class && !cov.bads.contains(cov.names[i]) )
        {
            picks[k] = i;
            ++k;
        }
    picks.conservativeResize(1, k);

    return picks;
}


/***********************************************************************************************
 *   SELECT APPROPRIATE CHANNELS
 **********************************************************************************************/

Eigen::MatrixXd selectAppropriateChannels_Matti(const Eigen::MatrixXd& data,
                                                const QList<FIFFLIB::FiffChInfo>& chs,
                                                QStringList picks)
{
    //chs = NULL;
    int j,k;
    Eigen::MatrixXd cov = Eigen::MatrixXd(picks.size(), picks.size());
//    double *cov_diag = NULL;
//    int   *is_meg = NULL;
    Eigen::MatrixXd res;

    //  ToDo: Take that error seriously
    if (picks.size() == 0)
    {
        std::cout << "No channels specified for picking in mne_pick_chs_cov_omit" << std::endl;
    }

    //  ToDo: Take that error seriously
    if (chs.size() == 0)
    {
        std::cout << "No names in covariance matrix. Cannot do picking." << std::endl;
    }

    //  setting up the pick vector, containing the indices, of the appropriate channels
    std::vector<int> pick = std::vector<int>(picks.size(), -1 );
    for (j = 0; j < picks.size(); ++j)
    {
        for (k = 0; k < chs.size(); ++k)
        {
            if ( chs[k].ch_name == picks[j] )
            {
                pick[j] = k;
                break;
            }
        }
    }

    //  check if picks was fully set up
    for (j = 0; j < picks.size(); ++j)
    {
        //  ToDo: Take that error seriously
        if (pick[j] < 0)
        {
          std::cout << "All desired channels not found in the covariance matrix (at least missing"
                    << picks[j].toStdString()
                    << ")." << std::endl;
        }
    }

//  Stuff in Mattis Code, that is irrelevant in this setting

    QStringList names;
    //  ToDo: Check this out
//    if (c->cov_diag) {
//    cov_diag = MALLOC(ncov,double);
//    for (j = 0; j < ncov; j++) {
//      cov_diag[j] = c->cov_diag[pick[j]];
//      names[j] = mne_strdup(c->names[pick[j]]);
//    }
//    }
//    else
    {
        for (j = 0; j < picks.size(); j++)
        {
            names.append(chs[pick[j]].ch_name);
            for (k = 0; k <= j; k++)
            {
                //  Stuff in Mattis Code, that cannot happen in this data structure
                cov(j,k) = data(pick[j],pick[k]);
                cov(k,j) = data(pick[j],pick[k]);
                //  Stuff in Mattis Code that is irrelevant in this setting
            }
        }
    }

    return cov;
//    res = mne_new_cov(c->kind,ncov,names,cov,cov_diag);

//    res->bads = mne_dup_name_list(c->bads,c->nbad);
//    res->nbad = c->nbad;
//    res->proj = mne_dup_proj_op(c->proj);
//    res->sss  = mne_dup_sss_data(c->sss);

//    if (c->ch_class) {
//    res->ch_class = MALLOC(res->ncov,int);
//    for (k = 0; k < res->ncov; k++)
//      res->ch_class[k] = c->ch_class[pick[k]];
//    }
//    FREE(pick);
//    FREE(is_meg);
//    return res;
}


Eigen::MatrixXd selectAppropriateChannels_Thomas(const Eigen::MatrixXd& data, Eigen::RowVectorXi picks)
{
    Eigen::MatrixXd pickedData = Eigen::MatrixXd(picks.size(), data.cols());
    for( int i = 0; i != picks.size(); ++i )
    {
        for( int j = 0; j != data.cols(); ++j )
        {
            pickedData(i,j) = data(picks[i],j);
        }
    }

    return pickedData;
}

/********************************************************************************************************
 * Function - cumpute_ssp_vectors
 * Here everythin is added together
 *******************************************************************************************************/
 void  compute_ssp_vectors( FIFFLIB::FiffCov& cov,
                            QString tag,  //  ToDo: check if necessary, supposed to be used for labelling
                            FIFFLIB::fiff_int_t ch_class,
                            FIFFLIB::fiff_int_t ncomp,
                            QStringList browse_names    // ToDo: Delete this, when a working solution was found
                            )
{
    //  the matrices cov_Matti and cov_Thomas are just for testing the two different strategies against each other
     FIFFLIB::FiffCov cov_Matti = cov;
     FIFFLIB::FiffCov cov_Thomas = cov;
    //  ToDo: maybe some additional variables needed

    //  ToDo Check that cov.ch_class contains the necessary data,
    //  esp. cov.ch_class != empty


    /**********************************************************************************************************************
     *  CREATE CHANNEL SELECTION
     *  Check each channel, if its kind == 'ch_class'  //   ToDo: Check if kind is right here
     *  and if its name is in the bads
     **********************************************************************************************************************/
    std::cout << "Setting up pick list..." << std::flush ;

    QStringList picks_Matti = createChannelSelection_Matti( cov_Matti, ch_class);
    Eigen::RowVectorXi picks_Thomas = createChannelSelection_Thomas( cov_Thomas, ch_class);
    //  Debugging: test if both versions yield equal results
    {
        FIFFLIB::fiff_int_t difference = 0;
        if(picks_Matti.size() != picks_Thomas.size() )
        {
            std::cout << "ERROR at creation of channel selection! - Different Dimensionality" << std::endl;
        }
        for( int i = 0; i!= picks_Matti.size(); ++i )
        {
            if( cov.names[picks_Thomas[i]] != picks_Matti[i] )
            {
                ++difference;
                std::cout << "Index " << i << ": "
                          << cov.names[picks_Thomas[i]].toStdString()
                          << " vs. " << picks_Matti[i].toStdString() << std::endl;
            }
        }

        if( difference != 0 )
        {
            std::cout << "ERROR at creation of channel selection! - Different channels picked" << std::endl;
        }

        //  ToDo: Delete this, when a working solution was found
        {
        //  check if same channels as in computed data ware picked
        for( int i = 0; i!= picks_Matti.size(); ++i )
        {
            if( !browse_names.contains(picks_Matti[i]) )
            {
                std::cout << picks_Matti[i].toStdString() << " was not picked in mne_browse_raw!"
                          << std::endl;
            }
        }
        }
     }

    std::cout << "DONE - dimension = " << picks_Matti.size() << std::endl;

    //  ToDo: Check, if picks is empty, then return the empty Projs

    //  SELECT APPROPRIATE CHANNELS

    //Eigen::MatrixXd pickedData_Matti = selectAppropriateChannels_Matti(data, chs, picks_Matti);
    Eigen::MatrixXd pickedData_Thomas = selectAppropriateChannels_Thomas(cov.data, picks_Thomas);

    //  COMPUTE EIGENVALUE DECOMPOSITION

    //  Put Eigenvectors into a FIFF::Proj
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
     * Read the empty room file into data
     **********************************************************************************************************/

    Eigen::RowVectorXi picks;
    QList<FIFFLIB::FiffProj> projectionVectors;
    FIFFLIB::FiffCov covarianceMatrix;
    //  ToDo: Delete projNames, when a working solution was found
    QStringList projNames;

    // the {} are used to delete unnecessary data from the stack after this
    {
        //  location of the file that will be read
        QFile t_fileEmptyRoom("./MNE-sample-data/MEG/test_data_ssp/151015_131148_4884471_Empty_room_raw.fif");
        //QFile t_fileEmptyRoom("./MNE-sample-data/MEG/test_data_ssp/Empty_Room_new_raw.fif");

        bool in_samples = true;

        //   Setup for reading the raw data
        FIFFLIB::FiffRawData raw(t_fileEmptyRoom);
        FIFFLIB::fiff_int_t from = 0;
        FIFFLIB::fiff_int_t to = 1000;


        //
        //  Set up pick list will be omitted here.
        //  Everything is picked and the selection will be done in the computeSspVectors-function
        picks = Eigen::RowVectorXi(raw.info.ch_names.size() );
        for( int i = 0; i != picks.cols(); ++i )
        {
            picks[i] = 1;
        }

        //  delete projections of raw file
        //  otherwise the projection will be carried out while reading ot the file
        raw.proj.resize(0,0);

        //   Read a data segment
        //   times output argument is optional
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

    /******************************************************************************************************************************
     * prepare covarianceMatrix, so far by hand
     *****************************************************************************************************************************/

        covarianceMatrix.kind = -1;             //  ToDp: This is not thought through, this is just for testing
        covarianceMatrix.diag = false;          //  ToDp: This is not thought through, this is just for testing
        covarianceMatrix.dim = data.rows();         /**< Dimension of the covariance (dim x dim). */
        covarianceMatrix.names = raw.info.ch_names;          /**< Channel names. */
        covarianceMatrix.data = data * data.transpose();          /**< Covariance data */
        covarianceMatrix.projs.clear();              /**< List of available ssp projectors. */
        covarianceMatrix.bads = raw.info.bads;       /**< List of bad channels. */
        //fiff_int_t nfree;         //  ToDo: Check what this is needed for and what it is supposed to be
                                    //  It probably is covMat.dim - #{eig != 0}, or the exact opposite!
                                    //  AKS SOMEBODY WHO KNOWS!
        //VectorXd eig;             //  This is calculated in computeSsp...
        //MatrixXd eigvec;          //  This is calculated in computeSsp...
        //  Set up covarianceMatrix.ch_class
        for ( int k = 0; k != raw.info.ch_names.size(); ++k )
        {
            covarianceMatrix.ch_class.push_back(raw.info.chs[k].kind);
        }

        //ToDo: Delete this, when a working solution was found
        {
        for( int i = 0; i != raw.info.projs[0].data->col_names.size(); ++i )
        {
            projNames.push_back(raw.info.projs[0].data->col_names[i]);
        }
        }

        t_fileEmptyRoom.close();
    }

    /*************************************************************************************************************************************************
     *  Compute the SSP vectors, using the covariance Matrix - almost all necessary data are in there
     * ToDo: Check, why dimensions don't fit with the already computed values - our result yields 1 channel more
     * The reason could very possibly lie in the selection by chs[k].kind == ch_class
     ************************************************************************************************************************************************/
    //  ToDo: Ask Lorenz about channel 'MEG348' - it's not picked in mne_browse raw, but I see no reason for that...
    compute_ssp_vectors( covarianceMatrix, "This is just a Test" , 1, 20, projNames /* ToDo:Delete projNames, when a working solution was found*/ );

    int temporaryStopper;
    std::cout << "Insert any number, to continue" << std::endl;
    std::cin >> temporaryStopper;

    /*************************************************************************************************************************************************
     *  Writing the projection vectors into a file with spontaneous data
     ************************************************************************************************************************************************/

    //  locations of the file to read and where to write the new file
    //  QFile t_fileIn("./MNE-sample-data/MEG/test_data_ssp/151015_151137_4884471_Spontaneous_raw.fif");
    QFile t_fileIn("./MNE-sample-data/MEG/test_data_ssp/151015_131148_4884471_Empty_room_raw.fif");
    //QFile t_fileIn("./MNE-sample-data/MEG/test_data_ssp/Empty_Room_new_raw.fif");

    QFile t_fileOut("./MNE-sample-data/MEG/test_data_ssp/ssp_output_spontaneous_with_new_projs_raw.fif");

    //
    //   Setup for reading the raw data
    //


    FIFFLIB::FiffRawData rawIn(t_fileIn);

    //use old pick list!

    /**************************************************************************************************************
     * Check if the new SSP vectors are different from Mattis vectors
     *************************************************************************************************************/

    int count = 0;
    if ( rawIn.info.projs.size() == projectionVectors.size() )
    {
        for( FIFFLIB::fiff_int_t i = 0; i != 2/*i != projectionVectors.size()*/; ++i )
        {
            Eigen::MatrixXd vecMatti = rawIn.info.projs[i].data->data;

            for( FIFFLIB::fiff_int_t j = 0; j != projectionVectors.size(); ++j )
            {
                Eigen::MatrixXd vecThomas = projectionVectors[j].data->data;
                Eigen::MatrixXd diff = vecMatti;
                diff -= vecThomas;
                std::cout << "Norm of difference VM(" << i << ") -VT(" << j << ") = " << diff.norm()
                          //<< "\nScalarproduct of the two vectors = " << vecMatti * vecThomas.transpose()
                          << std::endl;
            }
//            Eigen::MatrixXd vecThomas = projectionVectors[i].data->data;
//            if (vecMatti.cols() != vecThomas.cols() || vecMatti.rows()  != vecThomas.rows() )
//            {
//                std::cout << "WARNING: different sizes of vectors:/n"
//                          << "Mattis vector: " << vecMatti.rows() << "x" << vecMatti.cols() << "/n"
//                          << "Thomas vector: " << vecThomas.rows() << "x" << vecThomas.cols() << std::endl;
//                for( int k = 0; k != projectionVectors[i].data->col_names.size(); ++k )
//                {
//                    if( rawIn.info.projs[i].data->col_names.contains(projectionVectors[i].data->col_names[k]) )
//                        ++count;
//                    else
//                    {
//                        std::cout << "Not found: " << projectionVectors[i].data->col_names[k].toStdString() << std::endl;
//                    }
//                    if( rawIn.info.projs[i].data->col_names[k] != projectionVectors[i].data->col_names[k])
//                    {
//                        std::cout <<  "WARNING: Differenc Channel names in equal components:"
//                                  <<  rawIn.info.projs[i].data->col_names[k].toStdString()
//                                  <<  " vs. "
//                                  <<  projectionVectors[i].data->col_names[k].toStdString()
//                                  << std::endl;
//                    }
//                }
//            }
//            else
//            {
//                std::cout << "Norm of Mattis vector = " << vecMatti.norm()
//                          << ", Norm of Thomas vector = " << vecThomas.norm()
//                          << std::endl;
//                Eigen::MatrixXd diff = vecMatti;
//                diff -= vecThomas;
//                std::cout << "Norm of difference = " << diff.norm()
//                          << "\nScalarproduct of the two vectors = " << vecMatti * vecThomas.transpose()
//                          << std::endl;
//            }
        }
    }
    else
        std::cout << "WARNING: different number of projections" << std::endl;



    /*********************************************************************************************************
     * continue with the writing
     ********************************************************************************************************/

    rawIn.proj.resize(0,0);

    rawIn.info.projs = projectionVectors;
    Eigen::RowVectorXd cals;

    Eigen::MatrixXd data;
    Eigen::MatrixXd times;


    FIFFLIB::FiffStream::SPtr outfid = FIFFLIB::Fiff::start_writing_raw(t_fileOut,rawIn.info, cals/*, picks*/);

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
















