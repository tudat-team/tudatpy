/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <stdexcept>

#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/Astrodynamics/EarthOrientation/readAmplitudeAndArgumentMultipliers.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to read the amplitudes and fundamental argument multiplers for tidal corrections
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > readAmplitudesAndFundamentalArgumentMultipliers(
        const std::string amplitudesFile,
        const std::string fundamentalArgumentMultipliersFile,
        const double minimumAmplitude )
{

    // Read amplitudes and fundamental argument multipliers into matrices.
    Eigen::MatrixXd amplitudesRaw = input_output::readMatrixFromFile( amplitudesFile );
    Eigen::MatrixXd fundamentalArgumentMultipliersRaw = input_output::readMatrixFromFile( fundamentalArgumentMultipliersFile );

    // Check whether amplitudes and fundamental argument multipliers matrices have same number of rows
    if( amplitudesRaw.rows( ) != fundamentalArgumentMultipliersRaw.rows( ) )
    {
        throw std::runtime_error( "Amplitude and argument multipler files contain unequal set of entries" );
    }

    // Check whether number of columns in fundamental argument multiplier matrix is equal to 6.
    if( fundamentalArgumentMultipliersRaw.cols( ) != 6 )
    {
        throw std::runtime_error( "Number of columns in fundamental argument multipler matrix not equal to 6" );
    }


    // If acceptance amplitude is larger than 0.0 check RSS of amplitudes for each entry.
    Eigen::MatrixXd fundamentalArgumentMultipliers, amplitudes;
    if( minimumAmplitude > 0.0 )
    {
        // Filter raw data and remove entries with RSS amplitude that is too low
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > filteredData =
                filterRawDataForAmplitudes( fundamentalArgumentMultipliersRaw, amplitudesRaw, minimumAmplitude );
        amplitudes = filteredData.first;
        fundamentalArgumentMultipliers = filteredData.second;
    }
    // If all amplitudes are to be expected, raw data is same as return data.
    else
    {
        amplitudes = amplitudesRaw;
        fundamentalArgumentMultipliers = fundamentalArgumentMultipliersRaw;
    }

    // return amplitudes and fundamental argument multipliers.
    return std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( amplitudes, fundamentalArgumentMultipliers );
}

//! Function to filter tidal corrections based on RSS amplitude criteria.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > filterRawDataForAmplitudes(
        const Eigen::MatrixXd rawAmplitudes,
        const Eigen::MatrixXd rawFundamentalArgumentMultipliers,
        const double minimumAmplitude )
{
    // Get number of amplitudes per entry.
    int numberOfAmplitudes = rawAmplitudes.cols( );

    // Initialize return matrices
    int numberOfAcceptedAmplitudes = 0;
    Eigen::MatrixXd filteredAmplitudes;
    filteredAmplitudes.resize( numberOfAcceptedAmplitudes, numberOfAmplitudes );
    Eigen::MatrixXd filteredFundamentalArgumentMultipliers;
    filteredFundamentalArgumentMultipliers.resize( numberOfAcceptedAmplitudes, 6 );

    // Iterate over all entries, calculate RSS amplitude and accept or reject entry.
    double rssAmplitude = 0.0;
    for( int i = 0; i < rawAmplitudes.rows( ); i++ )
    {
        // Calculate RSS of amplitude.
        rssAmplitude = 0.0;
        for( int j =0; j < numberOfAmplitudes; j++ )
        {
            rssAmplitude += rawAmplitudes( i, j ) * rawAmplitudes( i, j );
        }
        rssAmplitude = std::sqrt( rssAmplitude );

        // If RSS amplitude is sufficiently large, accept value.
        if( rssAmplitude > minimumAmplitude )
        {
            if( numberOfAcceptedAmplitudes == 0 )
            {
                numberOfAcceptedAmplitudes++;
                filteredAmplitudes.resize( 1, numberOfAmplitudes );
                filteredAmplitudes = rawAmplitudes.block( i, 0, 1, numberOfAmplitudes );

                filteredFundamentalArgumentMultipliers.resize( 1, 6 );
                filteredFundamentalArgumentMultipliers = rawFundamentalArgumentMultipliers.block( i, 0, 1, 6 );
            }
            else
            {
                numberOfAcceptedAmplitudes++;
                filteredAmplitudes.conservativeResize( numberOfAcceptedAmplitudes, numberOfAmplitudes );
                filteredAmplitudes.block( numberOfAcceptedAmplitudes - 1, 0, 1, numberOfAmplitudes ) =
                        rawAmplitudes.block( i, 0, 1, numberOfAmplitudes );

                filteredFundamentalArgumentMultipliers.conservativeResize( numberOfAcceptedAmplitudes, 6 );
                filteredFundamentalArgumentMultipliers.block( numberOfAcceptedAmplitudes - 1, 0, 1, 6 )  =
                        rawFundamentalArgumentMultipliers.block( i, 0, 1, 6 );
            }
        }
    }

    return std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( filteredAmplitudes, filteredFundamentalArgumentMultipliers );
}

}

}

