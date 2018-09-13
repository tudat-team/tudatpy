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

#define BOOST_TEST_MAIN

#include <map>
#include <limits>
#include <iostream>

#include <Eigen/Core>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_basic_statistics )

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testSampleMean )
{
    // Test computation of sample mean on finite population using unbiased estimators.
    // The expected values are computed using the Microsoft Excel the AVERAGE( ) function.

    // Declare vector of sample data.
    std::vector< double > sampleData;

    // Populate vector with sample data.
    sampleData.push_back( 2.5 );
    sampleData.push_back( 6.4 );
    sampleData.push_back( 8.9 );
    sampleData.push_back( 12.7 );
    sampleData.push_back( 15.0 );

    // Set expected sample mean.
    double expectedSampleMean = 9.1;

    // Compute sample mean.
    double computedSampleMean = statistics::computeSampleMean( sampleData );

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSampleMean, expectedSampleMean,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample variance is computed correctly.
BOOST_AUTO_TEST_CASE( testSampleVariance )
{
    // Test computation of sample variance on finite population using unbiased estimators.
    // The expected values are computed using the Microsoft Excel the VAR( ) function.

    // Declare vector of sample data.
    std::vector< double > sampleData;

    // Populate vector with sample data.
    sampleData.push_back( 2.5 );
    sampleData.push_back( 6.4 );
    sampleData.push_back( 8.9 );
    sampleData.push_back( 12.7 );
    sampleData.push_back( 15.0 );

    // Declare expected sample variance.
    double expectedSampleVariance = 24.665;

    // Compute sample variance.
    double computedSampleVariance = statistics::computeSampleVariance( sampleData );

    // Check if computed sample variance matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSampleVariance, expectedSampleVariance,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if moving average is computed correctly. Results compared with MATLAB movmean function.
BOOST_AUTO_TEST_CASE( testMovingAverage )
{
    // Test case with simple VectorXd
    {
        // Populate vector with randomly distributed data points
        Eigen::VectorXd vectorOfPoints;
        vectorOfPoints.resize( 15 );
        vectorOfPoints << 0.93329860101525, 3.93372816267124, 1.35032100135611, 2.97099423629127, 3.18245216750598, 1.43494398584927,
                0.915460520182276, 2.60394635060288, 2.09834777464011, 3.04137361348961, 0.265830887303261, 2.96918626998768,
                1.23234701262448, 2.42638755740894, 1.6271912582765;

        // Compute expected moving average with MATLAB
        Eigen::VectorXd expectedMovingAverage;
        expectedMovingAverage.resize( 15 );
        expectedMovingAverage << 2.07244925501420, 2.29708550033347, 2.47415883376797, 2.57448791073478, 1.97083438223698, 2.22155945208634,
                2.04703015975610, 2.01881444895283, 1.78499182924363, 2.19573697920471, 1.92141711160903, 1.98702506816280, 1.70418859712017,
                2.06377802457440, 1.76197527610331;

        // Compute moving average with built-in function
        Eigen::VectorXd computedMovingAverage = statistics::computeMovingAverage( vectorOfPoints, 5 );

        // Compare results
        double tolerance = 1.0e2 * std::numeric_limits< double >::epsilon( );
        for ( unsigned int i = 0; i < 15; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( computedMovingAverage[ i ], expectedMovingAverage[ i ],
                                        tolerance );
        }
    }

    // Test case with vector of data
    {
        // Populate map with randomly distributed data points
        std::vector< Eigen::Vector3d > sampleData;
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 1.65235588866137, 2.32705996717709, 2.08263350423676 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 4.00607711081905, 2.34909226340225, 3.25705615743397 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 0.05562219359578, 0.67821147860744, 2.92482593349371 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 3.00004984907525, 1.94508108539059, 3.91112726565386 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 2.59458369740905, 2.35020117387454, 4.25025122830500 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 1.92978945855772, 1.23976325705858, 0.30963889688877 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 1.34844635824972, 4.19210187053127, 1.38816961132219 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 0.97553806336408, 0.05115282310110, 3.02049801445265 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << 2.86171630239342, 2.00116208348351, 1.92916278683952 ).finished( ) );
        sampleData.push_back( ( Eigen::VectorXd( 3 ) << -1.48628392070328, 2.58117232267592, -1.19243491996590 ).finished( ) );

        // Compute expected moving average with MATLAB
        std::map< double, Eigen::Vector3d > expectedMovingAverage;
        expectedMovingAverage[ 0.0 ] = ( Eigen::VectorXd( 3 ) << 1.90468506435874, 1.78478790306226, 2.75483853172148 ).finished( );
        expectedMovingAverage[ 1.0 ] = ( Eigen::VectorXd( 3 ) << 2.17852626053786, 1.82486119864434, 3.04391071520457 ).finished( );
        expectedMovingAverage[ 2.0 ] = ( Eigen::VectorXd( 3 ) << 2.26173774791210, 1.92992919369038, 3.28517881782466 ).finished( );
        expectedMovingAverage[ 3.0 ] = ( Eigen::VectorXd( 3 ) << 2.31722446189137, 1.71246985166668, 2.93057989635506 ).finished( );
        expectedMovingAverage[ 4.0 ] = ( Eigen::VectorXd( 3 ) << 1.78569831137750, 2.08107177309248, 2.55680258713270 ).finished( );
        expectedMovingAverage[ 5.0 ] = ( Eigen::VectorXd( 3 ) << 1.96968148533116, 1.95566004199122, 2.57593700332449 ).finished( );
        expectedMovingAverage[ 6.0 ] = ( Eigen::VectorXd( 3 ) << 1.94201477599480, 1.96687624160980, 2.17954410756163 ).finished( );
        expectedMovingAverage[ 7.0 ] = ( Eigen::VectorXd( 3 ) << 1.12584125237233, 2.01307047137008, 1.09100687790745 ).finished( );
        expectedMovingAverage[ 8.0 ] = ( Eigen::VectorXd( 3 ) << 0.92485420082599, 2.20639727494795, 1.28634887316211 ).finished( );
        expectedMovingAverage[ 9.0 ] = ( Eigen::VectorXd( 3 ) << 0.78365681501807, 1.54449574308685, 1.25240862710876 ).finished( );

        // Compute moving average with built-in function
        std::vector< Eigen::Vector3d > computedMovingAverage = statistics::computeMovingAverage( sampleData, 5 );

        // Compare results
        unsigned int i = 0;
        double tolerance = 1.0e2 * std::numeric_limits< double >::epsilon( );
        for ( std::map< double, Eigen::Vector3d >::const_iterator mapIterator = expectedMovingAverage.begin( );
              mapIterator != expectedMovingAverage.end( ); mapIterator++, i++ )
        {
            for ( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( computedMovingAverage.at( i )[ j ], mapIterator->second[ j ],
                                            tolerance );
            }
        }
    }

    // Test case with map of data
    {
        // Populate map with randomly distributed data points
        std::map< double, Eigen::VectorXd > sampleData;
        sampleData[ 0.0 ] = ( Eigen::VectorXd( 3 ) << 1.65235588866137, 2.32705996717709, 2.08263350423676 ).finished( );
        sampleData[ 1.0 ] = ( Eigen::VectorXd( 3 ) << 4.00607711081905, 2.34909226340225, 3.25705615743397 ).finished( );
        sampleData[ 2.0 ] = ( Eigen::VectorXd( 3 ) << 0.05562219359578, 0.67821147860744, 2.92482593349371 ).finished( );
        sampleData[ 3.0 ] = ( Eigen::VectorXd( 3 ) << 3.00004984907525, 1.94508108539059, 3.91112726565386 ).finished( );
        sampleData[ 4.0 ] = ( Eigen::VectorXd( 3 ) << 2.59458369740905, 2.35020117387454, 4.25025122830500 ).finished( );
        sampleData[ 5.0 ] = ( Eigen::VectorXd( 3 ) << 1.92978945855772, 1.23976325705858, 0.30963889688877 ).finished( );
        sampleData[ 6.0 ] = ( Eigen::VectorXd( 3 ) << 1.34844635824972, 4.19210187053127, 1.38816961132219 ).finished( );
        sampleData[ 7.0 ] = ( Eigen::VectorXd( 3 ) << 0.97553806336408, 0.05115282310110, 3.02049801445265 ).finished( );
        sampleData[ 8.0 ] = ( Eigen::VectorXd( 3 ) << 2.86171630239342, 2.00116208348351, 1.92916278683952 ).finished( );
        sampleData[ 9.0 ] = ( Eigen::VectorXd( 3 ) << -1.48628392070328, 2.58117232267592, -1.19243491996590 ).finished( );

        // Compute expected moving average with MATLAB
        std::map< double, Eigen::Vector3d > expectedMovingAverage;
        expectedMovingAverage[ 0.0 ] = ( Eigen::VectorXd( 3 ) << 1.90468506435874, 1.78478790306226, 2.75483853172148 ).finished( );
        expectedMovingAverage[ 1.0 ] = ( Eigen::VectorXd( 3 ) << 2.17852626053786, 1.82486119864434, 3.04391071520457 ).finished( );
        expectedMovingAverage[ 2.0 ] = ( Eigen::VectorXd( 3 ) << 2.26173774791210, 1.92992919369038, 3.28517881782466 ).finished( );
        expectedMovingAverage[ 3.0 ] = ( Eigen::VectorXd( 3 ) << 2.31722446189137, 1.71246985166668, 2.93057989635506 ).finished( );
        expectedMovingAverage[ 4.0 ] = ( Eigen::VectorXd( 3 ) << 1.78569831137750, 2.08107177309248, 2.55680258713270 ).finished( );
        expectedMovingAverage[ 5.0 ] = ( Eigen::VectorXd( 3 ) << 1.96968148533116, 1.95566004199122, 2.57593700332449 ).finished( );
        expectedMovingAverage[ 6.0 ] = ( Eigen::VectorXd( 3 ) << 1.94201477599480, 1.96687624160980, 2.17954410756163 ).finished( );
        expectedMovingAverage[ 7.0 ] = ( Eigen::VectorXd( 3 ) << 1.12584125237233, 2.01307047137008, 1.09100687790745 ).finished( );
        expectedMovingAverage[ 8.0 ] = ( Eigen::VectorXd( 3 ) << 0.92485420082599, 2.20639727494795, 1.28634887316211 ).finished( );
        expectedMovingAverage[ 9.0 ] = ( Eigen::VectorXd( 3 ) << 0.78365681501807, 1.54449574308685, 1.25240862710876 ).finished( );

        // Compute moving average with built-in function
        std::map< double, Eigen::VectorXd > computedMovingAverage = statistics::computeMovingAverage( sampleData, 5 );

        // Compare results
        double tolerance = 1.0e2 * std::numeric_limits< double >::epsilon( );
        for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = computedMovingAverage.begin( );
              mapIterator != computedMovingAverage.end( ); mapIterator++ )
        {
            for ( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( mapIterator->second[ i ], expectedMovingAverage[ mapIterator->first ][ i ],
                                            tolerance );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
