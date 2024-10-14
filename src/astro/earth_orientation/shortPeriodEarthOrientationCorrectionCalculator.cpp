/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to sum all the corrcetion terms.
template< >
double ShortPeriodEarthOrientationCorrectionCalculator< double >::sumCorrectionTerms( const Eigen::Vector6d& arguments )
{
    // Initialize corrections to zero.
    double currentCorrection = 0;
    double tideAngle = 0.0;

    // Iterate over all correction types
    for( unsigned int i = 0; i < argumentAmplitudes_.size( ); i++ )
    {
        Eigen::MatrixXd currentArgumentMultipliers = argumentMultipliers_.at( i );
        Eigen::MatrixXd currentAmplitudes = argumentAmplitudes_.at( i );

        // Iterte over all libration corrections.
        for( int j = 0; j < currentAmplitudes.rows( ); j++ )
        {
            // Calculate current phase angle of tide.
            tideAngle = ( arguments.transpose( ) *
                          currentArgumentMultipliers.block( j, 0, 1, 6 ).transpose( ) )( 0, 0 );

            currentCorrection += currentAmplitudes( j, 0 ) * std::sin( tideAngle ) +
                    currentAmplitudes( j, 1 ) * std::cos( tideAngle );
        }
    }


    return currentCorrection;
}

//! Function to sum all the corrcetion terms.
template< >
Eigen::Vector2d ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d >::sumCorrectionTerms(
        const Eigen::Vector6d& arguments )
{

    // Initialize corrections to zero.
    Eigen::Vector2d currentCorrection = Eigen::Vector2d::Zero( );
    double tideAngle = 0.0;

    // Iterate over all correction types
    for( unsigned int i = 0; i < argumentAmplitudes_.size( ); i++ )
    {
        Eigen::MatrixXd currentArgumentMultipliers = argumentMultipliers_.at( i );
        Eigen::MatrixXd currentAmplitudes = argumentAmplitudes_.at( i );

        // Iterte over all libration corrections.
        for( int j = 0; j < currentAmplitudes.rows( ); j++ )
        {
            // Calculate current phase angle of tide.
            tideAngle = ( arguments.transpose( ) *
                          currentArgumentMultipliers.block( j, 0, 1, 6 ).transpose( ) )( 0, 0 );

            // Calculate and add sine and cosine terms
            currentCorrection.x( ) += currentAmplitudes( j, 0 ) * std::sin( tideAngle ) +
                    currentAmplitudes( j, 1 ) * std::cos( tideAngle );
            currentCorrection.y( ) += currentAmplitudes( j, 2 ) * std::sin( tideAngle ) +
                    currentAmplitudes( j, 3 ) * std::cos( tideAngle );
        }
    }


    return currentCorrection;
}

//! Function to retrieve the default UT1 short-period correction calculator
std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > getDefaultUT1CorrectionCalculator(
        const double minimumAmplitude )
{

    return std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                1.0E-6, minimumAmplitude,
                std::vector< std::string >{
                    tudat::paths::getEarthOrientationDataFilesPath(  ) + "/utcLibrationAmplitudes.txt",
                    tudat::paths::getEarthOrientationDataFilesPath(  ) + "/utcOceanTidesAmplitudes.txt" },
                std::vector< std::string >{
                    tudat::paths::getEarthOrientationDataFilesPath(  ) + "/utcLibrationFundamentalArgumentMultipliers.txt",
                    tudat::paths::getEarthOrientationDataFilesPath(  ) + "/utcOceanTidesFundamentalArgumentMultipliers.txt" },
                std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ) );
}

//! Function to retrieve the default polar motion short-period correction calculator
std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > getDefaultPolarMotionCorrectionCalculator(
        const double minimumAmplitude )
{

    return std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                unit_conversions::convertArcSecondsToRadians< double >( 1.0E-6 ), minimumAmplitude,
                std::vector< std::string >{
                    tudat::paths::getEarthOrientationDataFilesPath(  ) +
                    "/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
                    tudat::paths::getEarthOrientationDataFilesPath(  ) +
                    "/polarMotionOceanTidesAmplitudes.txt", },
                std::vector< std::string >{
                    tudat::paths::getEarthOrientationDataFilesPath(  ) +
                    "/polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt",
                    tudat::paths::getEarthOrientationDataFilesPath(  ) +
                    "/polarMotionOceanTidesFundamentalArgumentMultipliers.txt" },
                std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ) );

}


}

}
