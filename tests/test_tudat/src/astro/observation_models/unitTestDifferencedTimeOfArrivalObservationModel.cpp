/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
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

BOOST_AUTO_TEST_SUITE( test_differenced_time_of_arrival )

BOOST_AUTO_TEST_CASE( testDifferencedTimeOfArrival )
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

    // Create observation settings
    std::shared_ptr< ObservationModelSettings > rangeObservableSettings1 = std::make_shared< ObservationModelSettings >(
            one_way_range,
            linkEnds1,
            lightTimeCorrectionSettings );

    std::shared_ptr< ObservationModelSettings > rangeObservableSettings2 = std::make_shared< ObservationModelSettings >(
            one_way_range,
            linkEnds2,
            lightTimeCorrectionSettings );

    for( int testTimeScale = 0; testTimeScale < 2; testTimeScale ++ )
    {
        std::shared_ptr< ObservationModelSettings > differencedObservableSettings = std::make_shared< DifferencedTimeOfArrivalObservationSettings >(
                differencedLinkEnds,
                lightTimeCorrectionSettings,
                ( testTimeScale == 0 ) ? basic_astrodynamics::tdb_scale : basic_astrodynamics::utc_scale );


        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, Time > > firstRangeObservationModel =
                ObservationModelCreator< 1, double, Time >::createObservationModel( rangeObservableSettings1, bodies );
        std::shared_ptr< ObservationModel< 1, double, Time > > secondRangeObservationModel =
                ObservationModelCreator< 1, double, Time >::createObservationModel( rangeObservableSettings2, bodies );
        std::shared_ptr< ObservationModel< 1, double, Time > > differencedTimeOfArrivalObservationModel =
                ObservationModelCreator< 1, double, Time >::createObservationModel( differencedObservableSettings, bodies );

        // Compute observation separately with two functions.
        Time receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
        std::vector< double > linkEndTimesDifferenced, linkEndTimesRange1, linkEndTimesRange2;
        std::vector< Eigen::Vector6d > linkEndStatesDifferenced, linkEndStatesRange1, linkEndStatesRange2;

        double differencedTimeOfArrival =
                differencedTimeOfArrivalObservationModel->computeObservationsWithLinkEndData(
                        receiverObservationTime, receiver, linkEndTimesDifferenced, linkEndStatesDifferenced )( 0 );


        double firstRange =
                firstRangeObservationModel->computeObservationsWithLinkEndData(
                        receiverObservationTime, receiver, linkEndTimesRange1, linkEndStatesRange1 )( 0 );


         double secondRange =
                 secondRangeObservationModel->computeObservationsWithLinkEndData(
                         linkEndTimesRange1.at( 0 ), transmitter, linkEndTimesRange2, linkEndStatesRange2 )( 0 );

        std::cout<<std::setprecision( 16 )<<linkEndTimesDifferenced.at( 0 )<<" "<<linkEndTimesDifferenced.at( 1 )<<" "<<linkEndTimesDifferenced.at( 2 )<<std::endl;
        std::cout<<linkEndTimesRange1.at( 0 )<<" "<<linkEndTimesRange1.at( 1 )<<std::endl;
        std::cout<<linkEndTimesRange2.at( 0 )<<" "<<linkEndTimesRange2.at( 1 )<<std::endl;

        BOOST_CHECK_CLOSE_FRACTION( linkEndTimesDifferenced.at( 0 ), linkEndTimesRange1.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        BOOST_CHECK_CLOSE_FRACTION( linkEndTimesDifferenced.at( 1 ), linkEndTimesRange1.at( 1 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStatesDifferenced.at( 0 ), linkEndStatesRange1.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStatesDifferenced.at( 1 ), linkEndStatesRange1.at( 1 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

        BOOST_CHECK_CLOSE_FRACTION( linkEndTimesDifferenced.at( 0 ), linkEndTimesRange2.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        BOOST_CHECK_CLOSE_FRACTION( linkEndTimesDifferenced.at( 2 ), linkEndTimesRange2.at( 1 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStatesDifferenced.at( 0 ), linkEndStatesRange2.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linkEndStatesDifferenced.at( 2 ), linkEndStatesRange2.at( 1 ), 2.0E-14 );

        std::cout<<( firstRange - secondRange ) / physical_constants::SPEED_OF_LIGHT<<" "<<differencedTimeOfArrival<<std::endl;
        std::cout<<( firstRange - secondRange ) / physical_constants::SPEED_OF_LIGHT - differencedTimeOfArrival<<std::endl;

        if( testTimeScale == 0 )
        {
            BOOST_CHECK_SMALL( ( firstRange - secondRange ) / physical_constants::SPEED_OF_LIGHT - differencedTimeOfArrival, 5.0E-13 );
        }
        else if( testTimeScale == 1 )
        {
            std::shared_ptr< TerrestrialTimeScaleConverter > defaultTimeConverter = createDefaultTimeConverter( );
            Time receptionTdbTime1 = receiverObservationTime;
            Time receptionTdbTime2 = receiverObservationTime - ( ( firstRange - secondRange ) / physical_constants::SPEED_OF_LIGHT );
            Time receptionUtcTime1 = defaultTimeConverter->getCurrentTime< Time >(
                    basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, receptionTdbTime1,
                    getApproximateDsnGroundStationPositions( ).at( "DSS-13" ) );
            Time receptionUtcTime2 = defaultTimeConverter->getCurrentTime< Time >(
                    basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, receptionTdbTime2,
                    getApproximateDsnGroundStationPositions( ).at( "DSS-43" ) );

            BOOST_CHECK_SMALL( static_cast< double >( receptionUtcTime1 - receptionUtcTime2 ) - differencedTimeOfArrival, 5.0E-13 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
