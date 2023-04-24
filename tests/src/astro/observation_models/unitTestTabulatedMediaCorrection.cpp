/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "tudat/basics/testMacros.h"

#include <numeric>

#include "tudat/astro/observation_models/corrections/tabulatedMediaCorrection.h"

namespace tudat
{
namespace unit_tests
{

using namespace observation_models;

BOOST_AUTO_TEST_SUITE( test_tabulated_media_corrections )

// Only testing wet simplified Chao model. Consistency of dry Chao model with Niell checked in following test
// Anyway, wet and dry Chao are identical but with different coefficients.
BOOST_AUTO_TEST_CASE( testSimplifiedChaoMappingFunction )
{

    // Estefan and Sovers (1994), A Comparative Survey of Current and Proposed Tropospheric Refraction-Delay Models for
    // DSN Radio Metric Data Calibration, JPL/NASA, 94-24
    // Appendix
    std::vector< double > chaoTabulatedWetMapping = {
            35.3955, 23.4935, 17.2402, 13.4500, 10.9967, 9.2827, 8.0235, 7.0621, 6.3054, 5.6951, 5.1929,
            4.7728, 4.4164, 4.1104, 3.8449, 3.6125, 3.4075, 3.2253, 3.0625, 2.9160, 2.7838, 2.6637, 2.5543
        };

    std::vector< double > chaoTabulatedMappingElevation( chaoTabulatedWetMapping.size( ) );
    for ( unsigned int i = 0; i < chaoTabulatedMappingElevation.size( ); ++i )
    {
        chaoTabulatedMappingElevation.at( i ) = ( i + 1.0 ) * mathematical_constants::PI / 180.0;
    }

    for ( unsigned int i = 0; i < chaoTabulatedWetMapping.size( ); ++i )
    {
        std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return chaoTabulatedMappingElevation.at( i ); };

        // isUplinkCorrection set to true. Irrelevant if true or false, since the elevation function returns a  fixed value
        SimplifiedChaoTroposphericMapping simplifiedChaoModel = SimplifiedChaoTroposphericMapping( elevationFunction, true );

        // Moyer (2000) says that error should be smaller than 1% for elevation larger than 1 deg. Not sure
        // why it needs higher tolerance below 3.5 deg
        double tolerance;
        if ( chaoTabulatedMappingElevation.at( i ) * 180.0 / mathematical_constants::PI > 3.5 )
        {
            tolerance = 0.01;
        }
        else
        {
            tolerance = 0.025;
        }
        // Arguments to computeWetTroposphericMapping are irrelevant
        BOOST_CHECK_CLOSE_FRACTION(
                chaoTabulatedWetMapping.at( i ),
                simplifiedChaoModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 ), tolerance );
    }

}

// Check consistency between Niell and Chao mapping
BOOST_AUTO_TEST_CASE( testNiellChaoMappingFunctionConsistency )
{

    std::vector< double > elevation( 17 );
    for ( unsigned int i = 0; i < elevation.size( ); ++i )
    {
        elevation.at( i ) = ( i + 2 ) * 5.0 * mathematical_constants::PI / 180.0;
    }

    // Geodetic position: [altitude, latitude, longitude]
    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction =
            [=]( double time ){ return ( Eigen::Vector3d( ) << 800.0, 35.0 * mathematical_constants::PI / 180.0, TUDAT_NAN ).finished( ); };

    // Comparing wet mapping
    for ( unsigned int i = 0; i < elevation.size( ); ++i )
    {
        std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return elevation.at( i ); };


        // isUplinkCorrection set to true. Irrelevant if true or false, since the elevation function returns a  fixed value
        SimplifiedChaoTroposphericMapping simplifiedChaoModel = SimplifiedChaoTroposphericMapping( elevationFunction, true );
        NiellTroposphericMapping niellModel = NiellTroposphericMapping(
                elevationFunction, groundStationGeodeticPositionFunction, true );

        double chaoMapping = simplifiedChaoModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );
        double niellMapping = niellModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );
        BOOST_CHECK_CLOSE_FRACTION( chaoMapping, niellMapping, 0.01 );
    }

    // Comparing dry mapping
    for ( unsigned int i = 0; i < elevation.size( ); ++i )
    {
        std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return elevation.at( i ); };

        // isUplinkCorrection set to true. Irrelevant if true or false, since the elevation function returns a  fixed value
        SimplifiedChaoTroposphericMapping simplifiedChaoModel = SimplifiedChaoTroposphericMapping( elevationFunction, true );
        NiellTroposphericMapping niellModel = NiellTroposphericMapping(
                elevationFunction, groundStationGeodeticPositionFunction, true );

        double chaoMapping = simplifiedChaoModel.computeDryTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );
        double niellMapping = niellModel.computeDryTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), 0.0, 0.0 );

        BOOST_CHECK_CLOSE_FRACTION( chaoMapping, niellMapping, 0.01 );
    }
}

// Compare values of Niell mapping function with
// Niell (1996), Global mapping functions for the atmosphere delay at radio wavelengths, Journal of Geophysics Research
BOOST_AUTO_TEST_CASE( testNiellMappingFunction )
{

    double elevation = 5.0 * mathematical_constants::PI / 180.0;
    double geodeticLatitude = 64.92 * mathematical_constants::PI / 180.0;
    double altitude = 132.0;

    double time1987 = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            1987, 1, 1, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    // Figure 3 of Niell (1996)
    std::vector< double > dryMapping = { 10.1875, 10.18, 10.1475, 10.1375, 10.185 };
    std::vector< double > dryTolerance = { 0.0025, 0.0025, 0.0025, 0.0025, 0.005 };
    std::vector< double > dryTime { time1987 + 0.1 * physical_constants::JULIAN_YEAR,
                                    time1987 + 0.2 * physical_constants::JULIAN_YEAR,
                                    time1987 + 0.4 * physical_constants::JULIAN_YEAR,
                                    time1987 + 0.7 * physical_constants::JULIAN_YEAR,
                                    time1987 + 1.0 * physical_constants::JULIAN_YEAR };

    // Figure 5 of Niell (1996)
    double wetMapping = 10.73;
    double wetTolerance = 0.01;
    std::vector< double > wetTime = dryTime;

    // Create mapping function
    std::function< double ( Eigen::Vector3d, double ) > elevationFunction =
                [=]( Eigen::Vector3d inertialVectorAwayFromStation, double time ){ return elevation; };
    // Geodetic position: [altitude, latitude, longitude]
    std::function< Eigen::Vector3d ( double ) > groundStationGeodeticPositionFunction =
            [=]( double time ){ return ( Eigen::Vector3d( ) << altitude, geodeticLatitude, TUDAT_NAN ).finished( ); };

    NiellTroposphericMapping niellModel = NiellTroposphericMapping(
            elevationFunction, groundStationGeodeticPositionFunction, true );

    // Test dry mapping
    for ( unsigned int i = 0; i < dryMapping.size( ); ++i )
    {
        double calculatedMapping = niellModel.computeDryTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), dryTime.at( i ), dryTime.at( i ) );

        BOOST_CHECK_SMALL( std::abs( calculatedMapping - dryMapping.at( i ) ), dryTolerance.at( i ) );
    }
    std::cout << std::endl;

    // Test wet mapping
    for ( unsigned int i = 0; i < wetTime.size( ); ++i )
    {
        double calculatedMapping = niellModel.computeWetTroposphericMapping(
                Eigen::Vector6d::Zero(), Eigen::Vector6d::Zero(), wetTime.at( i ), wetTime.at( i ) );

        BOOST_CHECK_SMALL( std::abs( calculatedMapping - wetMapping ), wetTolerance );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}