/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
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

#include "tudat/simulation/estimation_setup.h"

namespace tudat
{
namespace unit_tests
{

using namespace observation_models;

BOOST_AUTO_TEST_SUITE( test_solar_corona_corrections )

// Comparison with the correction from section 6.1 of:
// T. Morley and F. Budnik (2007) , EFFECTS ON SPACECRAFT RADIOMETRIC DATA AT SUPERIOR SOLAR CONJUNCTION,  20th
//      International Symposium on Space Flight Dynamics
BOOST_AUTO_TEST_CASE( testInversePowerSeriesCorrectionMorley )
{

    std::vector< double > sepAngle = { 10, 20, 30, 60, 90, 180 }; // deg
    std::vector< double > geocentricDistance = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 }; // AU

    // 2-way range correction
    std::vector< double > expectedRangeCorrection = {
            2.2, 19.0, 33.8, 35.8, 36.5, 36.8, // 10 deg
            2.0, 9.1, 14.6, 16.1, 16.7, 17.1, // 20 deg
            1.8, 5.8, 8.7, 9.8, 10.3, 10.6, // 30 deg
            1.3, 2.7, 3.5, 4.0, 4.3, 4.5, // 60 deg
            1.0, 1.7, 2.2, 2.5, 2.6, 2.8, // 90 deg
            0.7, 1.1, 1.3, 1.5, 1.6, 1.7 // 180 deg
    };

    double tolerance = 0.1;

    // Convert angles to radians
    for ( unsigned int i = 0; i < sepAngle.size( ); ++i )
        sepAngle.at( i ) *= mathematical_constants::PI / 180.0;

    // Convert geocentric distances to meter
    for ( unsigned int i = 0; i < geocentricDistance.size( ); ++i )
        geocentricDistance.at( i ) *= physical_constants::ASTRONOMICAL_UNIT;

    // Convert 2-way range correction to 1-way range (assuming it's identical both ways)
    for ( unsigned int i = 0; i < expectedRangeCorrection.size( ); ++i )
        expectedRangeCorrection.at( i ) *= 0.5;

    // Define frequency function. They say X-band, i.e. 7 to 8 GHz
    // Morley and Budnik (2007) say they use X-band
    // Wikipedia: X-band for uplink: 7.145 - 7.190 GHz, for downlink 8.4 - 8.45 GHz (https://en.wikipedia.org/wiki/Deep_space_bands)
    // Selected the frequency to make one of the points coincide
    std::function< double ( std::vector< FrequencyBands >, double ) > frequencyFunction =
            []( std::vector< FrequencyBands >, double ){ return 7.66e9; };

    // Set coefficients
    const std::vector< double > coefficients = { 1.3e14, 0.5e12 };
    const std::vector< double > integerPositiveExponents = { 6.0, 2.0 };
    const std::vector< double > doublePositiveExponents = { 6.0 + 1e-12, 2.0 + 1e-12 };

    // Define Sun state function
    std::function< Eigen::Vector6d ( double time ) > sunStateFunction =
            []( double time ){ return Eigen::Vector6d::Zero( ); };
    // Define Earth state
    Eigen::Vector6d earthState = Eigen::Vector6d::Zero( );
    earthState( 0 ) = 1.0 * physical_constants::ASTRONOMICAL_UNIT;

    InversePowerSeriesSolarCoronaCorrection coronaCorrectionAnalytical = InversePowerSeriesSolarCoronaCorrection(
            observation_models::n_way_range, sunStateFunction, frequencyFunction,
            coefficients, integerPositiveExponents );

    InversePowerSeriesSolarCoronaCorrection coronaCorrectionNumerical = InversePowerSeriesSolarCoronaCorrection(
            observation_models::n_way_range, sunStateFunction, frequencyFunction,
            coefficients, doublePositiveExponents );

    std::shared_ptr< ObservationAncilliarySimulationSettings > dummyAncillarySettings = std::make_shared<
            ObservationAncilliarySimulationSettings >( );
    dummyAncillarySettings->setAncilliaryDoubleVectorData( frequency_bands, { TUDAT_NAN } );

    for ( unsigned int i = 0; i < sepAngle.size( ); ++i )
    {
        for ( unsigned int j = 0; j < geocentricDistance.size( ); ++j )
        {
            int id = i * sepAngle.size( ) + j;

            Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
            spacecraftState( 0 ) = 1.0 * physical_constants::ASTRONOMICAL_UNIT - std::cos( sepAngle.at( i ) ) * geocentricDistance.at( j );
            spacecraftState( 1 ) = std::sin( sepAngle.at( i ) ) * geocentricDistance.at( j );

            double calculatedCorrectionAnalytical = coronaCorrectionAnalytical.calculateLightTimeCorrectionWithMultiLegLinkEndStates(
                    { earthState, spacecraftState }, { 0.0, 0.0 }, 0, dummyAncillarySettings ) * SPEED_OF_LIGHT;

            double calculatedCorrectionNumerical = coronaCorrectionNumerical.calculateLightTimeCorrectionWithMultiLegLinkEndStates(
                    { spacecraftState, earthState }, { 0.0, 0.0 }, 0, dummyAncillarySettings ) * SPEED_OF_LIGHT;

            BOOST_CHECK_SMALL( calculatedCorrectionAnalytical - expectedRangeCorrection.at( id ), tolerance );
            BOOST_CHECK_SMALL( calculatedCorrectionNumerical - expectedRangeCorrection.at( id ), tolerance );

        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}