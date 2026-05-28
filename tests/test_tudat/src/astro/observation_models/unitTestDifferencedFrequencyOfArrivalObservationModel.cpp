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

#include <cmath>
#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::earth_orientation;

BOOST_AUTO_TEST_SUITE( test_differenced_frequency_of_arrival )

BOOST_AUTO_TEST_CASE( testDifferencedFrequencyOfArrival )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    defaultBodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Set frequency of transmitter (Mars) to a time-varying X-Band ramp so UTC and TDB paths are distinguishable.
    double transmitterFrequency = 8.4e9;
    bodies.getBody( "Mars" )->getVehicleSystems( )->setTransmittedFrequencyCalculator(
            std::make_shared< ground_stations::PiecewiseLinearFrequencyInterpolator >(
                    std::vector< Time >{ initialEphemerisTime - buffer },
                    std::vector< Time >{ finalEphemerisTime + buffer },
                    std::vector< double >{ 1.0 },
                    std::vector< double >{ transmitterFrequency } ) );

    // Define link ends for observations.
    LinkEnds linkEnds1;
    linkEnds1[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    linkEnds1[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    LinkEnds linkEnds2;
    linkEnds2[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-43" );
    linkEnds2[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    LinkEnds differencedLinkEnds;
    differencedLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    differencedLinkEnds[ receiver2 ] = std::make_pair< std::string, std::string >( "Earth", "DSS-43" );
    differencedLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    double tdbDifferencedFrequencyOfArrival = std::numeric_limits< double >::quiet_NaN( );
    double tdbFirstFrequency = std::numeric_limits< double >::quiet_NaN( );
    double tdbSecondFrequency = std::numeric_limits< double >::quiet_NaN( );

    for( int testTimeScale = 0; testTimeScale < 2; testTimeScale++ )
    {
        std::shared_ptr< ObservationModelSettings > frequencyObservableSettings1 = oneWayDopplerMeasuredFrequencySettings(
                linkEnds1,
                lightTimeCorrectionSettings,
                ( testTimeScale == 0 ) ? basic_astrodynamics::tdb_scale : basic_astrodynamics::utc_scale );
        std::shared_ptr< ObservationModelSettings > frequencyObservableSettings2 = oneWayDopplerMeasuredFrequencySettings(
                linkEnds2,
                lightTimeCorrectionSettings,
                ( testTimeScale == 0 ) ? basic_astrodynamics::tdb_scale : basic_astrodynamics::utc_scale );
        std::shared_ptr< ObservationModelSettings > differencedObservableSettings = differencedFrequencyOfArrivalObservationSettings(
                differencedLinkEnds,
                lightTimeCorrectionSettings,
                ( testTimeScale == 0 ) ? basic_astrodynamics::tdb_scale : basic_astrodynamics::utc_scale,
                nullptr,
                std::make_shared< LightTimeConvergenceCriteria >( ) );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, Time > > firstFrequencyObservationModel =
                ObservationModelCreator< 1, double, Time >::createObservationModel( frequencyObservableSettings1, bodies );
        std::shared_ptr< ObservationModel< 1, double, Time > > secondFrequencyObservationModel =
                ObservationModelCreator< 1, double, Time >::createObservationModel( frequencyObservableSettings2, bodies );
        std::shared_ptr< ObservationModel< 1, double, Time > > differencedFrequencyOfArrivalObservationModel =
                ObservationModelCreator< 1, double, Time >::createObservationModel( differencedObservableSettings, bodies );

        // Compute observation separately with two functions.
        Time receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
        std::vector< double > linkEndTimesDifferenced, linkEndTimesRange1, linkEndTimesRange2;
        std::vector< Eigen::Vector6d > linkEndStatesDifferenced, linkEndStatesRange1, linkEndStatesRange2;

        double differencedFrequencyOfArrival = differencedFrequencyOfArrivalObservationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimesDifferenced, linkEndStatesDifferenced )( 0 );

        double firstFrequency = firstFrequencyObservationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimesRange1, linkEndStatesRange1 )( 0 );

        double secondFrequency = secondFrequencyObservationModel->computeObservationsWithLinkEndData(
                receiverObservationTime, receiver, linkEndTimesRange2, linkEndStatesRange2 )( 0 );

        std::cout << std::setprecision( 16 ) << linkEndTimesDifferenced.at( 0 ) << " " << linkEndTimesDifferenced.at( 1 ) << " "
                  << linkEndTimesDifferenced.at( 2 ) << " " << linkEndTimesDifferenced.at( 3 ) << std::endl;
        std::cout << linkEndTimesRange1.at( 0 ) << " " << linkEndTimesRange1.at( 1 ) << std::endl;
        std::cout << linkEndTimesRange2.at( 0 ) << " " << linkEndTimesRange2.at( 1 ) << std::endl;

        // Check first sub-model: transmitter and receiver times/states match
        BOOST_CHECK_CLOSE_FRACTION(
                linkEndTimesDifferenced.at( 0 ), linkEndTimesRange1.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        BOOST_CHECK_CLOSE_FRACTION(
                linkEndTimesDifferenced.at( 1 ), linkEndTimesRange1.at( 1 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStatesDifferenced.at( 0 ), linkEndStatesRange1.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStatesDifferenced.at( 1 ), linkEndStatesRange1.at( 1 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

        // Check second sub-model: transmitter and receiver times/states match
        BOOST_CHECK_CLOSE_FRACTION(
                linkEndTimesDifferenced.at( 2 ), linkEndTimesRange2.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        BOOST_CHECK_CLOSE_FRACTION(
                linkEndTimesDifferenced.at( 3 ), linkEndTimesRange2.at( 1 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                linkEndStatesDifferenced.at( 2 ), linkEndStatesRange2.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStatesDifferenced.at( 3 ), linkEndStatesRange2.at( 1 ), 2.0E-14 );

        std::cout << ( firstFrequency - secondFrequency ) << " " << differencedFrequencyOfArrival << std::endl;
        std::cout << ( firstFrequency - secondFrequency ) - differencedFrequencyOfArrival << std::endl;

        if( testTimeScale == 0 )
        {
            tdbDifferencedFrequencyOfArrival = differencedFrequencyOfArrival;
            tdbFirstFrequency = firstFrequency;
            tdbSecondFrequency = secondFrequency;
            BOOST_CHECK_SMALL( ( firstFrequency - secondFrequency ) - differencedFrequencyOfArrival, 5.0E-13 );
        }
        else if( testTimeScale == 1 )
        {
            std::shared_ptr< TerrestrialTimeScaleConverter > defaultTimeConverter = createDefaultTimeConverter( );
            Time receptionUtcTime = defaultTimeConverter->getCurrentTime< Time >(
                    basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, receiverObservationTime,
                    getApproximateDsnGroundStationPositions( ).at( "DSS-13" ) );
            Time receptionTdbTime = defaultTimeConverter->getCurrentTime< Time >(
                    basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, receptionUtcTime,
                    getApproximateDsnGroundStationPositions( ).at( "DSS-13" ) );

            BOOST_CHECK_SMALL( static_cast< double >( receptionTdbTime - receiverObservationTime ), 5.0E-13 );
            BOOST_CHECK_SMALL( ( firstFrequency - secondFrequency ) - differencedFrequencyOfArrival, 5.0E-13 );
            BOOST_CHECK_GT( std::fabs( differencedFrequencyOfArrival - tdbDifferencedFrequencyOfArrival ), 1.0E-7 );
            BOOST_CHECK_GT( std::fabs( firstFrequency - tdbFirstFrequency ), 1.0E-7 );
            BOOST_CHECK_GT( std::fabs( secondFrequency - tdbSecondFrequency ), 1.0E-7 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
