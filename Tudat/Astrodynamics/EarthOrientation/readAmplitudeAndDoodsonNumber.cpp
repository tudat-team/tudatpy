#include <iostream>
#include <boost/exception/all.hpp>

#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Astrodynamics/EarthOrientation/readAmplitudeAndDoodsonNumber.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to read the amplitudes and doodson multiplers for tidal corrections
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > readAmplitudesAndDoodsonMultipliers(
        const std::string amplitudesFile,
        const std::string doodsonMultipliersFile,
        const double minimumAmplitude )
{
    // Read amplitudes and doodson numbers into matrices.
    Eigen::MatrixXd amplitudesRaw = input_output::readMatrixFromFile( amplitudesFile );
    Eigen::MatrixXd doodsonMultipliersRaw = input_output::readMatrixFromFile( doodsonMultipliersFile );

    // Check whether amplitudes and Doodson multipliers matrices have same number of rows
    if( amplitudesRaw.rows( ) != doodsonMultipliersRaw.rows( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Amplitude and doodson multipler files contain unequal set of entries" ) ) );
    }

    // Check whether number of columns in Doodson multiplier matrix is equal to 6.
    if( doodsonMultipliersRaw.cols( ) != 6 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Number of columns in Doodson multipler matrix not equal to 6" ) ) );
    }


    // If acceptance amplitude is larger than 0.0 check RSS of amplitudes for each entry.
    Eigen::MatrixXd doodsonMultipliers, amplitudes;
    if( minimumAmplitude > 0.0 )
    {
        // Filter raw data and remove entries with RSS amplitude that is too low
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > filteredData =
                filterRawDataForAmplitudes( doodsonMultipliersRaw, amplitudesRaw, minimumAmplitude );
        amplitudes = filteredData.first;
        doodsonMultipliers = filteredData.second;
    }
    // If all amplitudes are to be expected, raw data is same as return data.
    else
    {
        amplitudes = amplitudesRaw;
        doodsonMultipliers = doodsonMultipliersRaw;
    }

    // return amplitudes and Doodson multipliers.
    return std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( amplitudes, doodsonMultipliers );
}

//! Function to filter tidal corrections based on RSS amplitude criteria.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > filterRawDataForAmplitudes(
        const Eigen::MatrixXd rawAmplitudes,
        const Eigen::MatrixXd rawDoodsonMultipliers,
        const double minimumAmplitude )
{
    // Get number of amplitudes per entry.
    int numberOfAmplitudes = rawAmplitudes.cols( );

    // Initialize return matrices
    int numberOfAcceptedAmplitudes = 0;
    Eigen::MatrixXd filteredAmplitudes;
    filteredAmplitudes.resize( numberOfAcceptedAmplitudes, numberOfAmplitudes );
    Eigen::MatrixXd filteredDoodsonMultipliers;
    filteredDoodsonMultipliers.resize( numberOfAcceptedAmplitudes, 6 );

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

                filteredDoodsonMultipliers.resize( 1, 6 );
                filteredDoodsonMultipliers = rawDoodsonMultipliers.block( i, 0, 1, 6 );
            }
            else
            {
                numberOfAcceptedAmplitudes++;
                filteredAmplitudes.conservativeResize( numberOfAcceptedAmplitudes, numberOfAmplitudes );
                filteredAmplitudes.block( numberOfAcceptedAmplitudes - 1, 0, 1, numberOfAmplitudes ) =
                        rawAmplitudes.block( i, 0, 1, numberOfAmplitudes );

                filteredDoodsonMultipliers.conservativeResize( numberOfAcceptedAmplitudes, 6 );
                filteredDoodsonMultipliers.block( numberOfAcceptedAmplitudes - 1, 0, 1, 6 )  =
                        rawDoodsonMultipliers.block( i, 0, 1, 6 );
            }
        }
    }

    return std::pair< Eigen::MatrixXd, Eigen::MatrixXd >( filteredAmplitudes, filteredDoodsonMultipliers );
}

}

}

