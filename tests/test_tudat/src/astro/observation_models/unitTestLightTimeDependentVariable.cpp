/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;

BOOST_AUTO_TEST_SUITE( test_light_time_dependent_variable )

SystemOfBodies createTestBodies( )
{
    std::vector< std::string > bodyNames{ "Earth", "Moon", "Sun" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth" );

    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings =
            keplerEphemerisSettings( keplerElements, 0.0, getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    return bodies;
}

std::shared_ptr< ObservationCollection<> > simulateWithDependentVariables(
        const ObservableType observableType,
        const SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const std::vector< double >& observationTimes,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariablesList,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& correctionSettings =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) )
{
    std::shared_ptr< ObservationModelSettings > observationSettings =
            std::make_shared< ObservationModelSettings >( observableType, linkEnds, correctionSettings );
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
            createObservationSimulators( std::vector< std::shared_ptr< ObservationModelSettings > >{ observationSettings }, bodies );

    std::shared_ptr< TabulatedObservationSimulationSettings<> > simulationSettings =
            std::make_shared< TabulatedObservationSimulationSettings<> >( observableType,
                                                                          linkEnds,
                                                                          observationTimes,
                                                                          receiver,
                                                                          std::vector< std::shared_ptr< ObservationViabilitySettings > >{},
                                                                          nullptr,
                                                                          ancillarySettings );
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput{ simulationSettings };
    addDependentVariablesToObservationSimulationSettings( measurementSimulationInput, dependentVariablesList, bodies );

    return simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );
}

double getScalarDependentVariable( const std::shared_ptr< ObservationCollection<> >& collection,
                                   const std::shared_ptr< ObservationDependentVariableSettings >& settings,
                                   const double time )
{
    std::map< double, Eigen::VectorXd > history = collection->getDependentVariableHistory( settings );
    BOOST_REQUIRE( history.count( time ) == 1 );
    BOOST_REQUIRE_EQUAL( history.at( time ).size( ), 1 );
    return history.at( time )( 0 );
}

//! One-way light time is read from the solved leg cache, rather than reconstructed by subtracting
//! two large absolute epochs. A late epoch makes loss of precision in the latter approach visible.
BOOST_AUTO_TEST_CASE( testOneWayLightTimeDependentVariable )
{
    loadStandardSpiceKernels( );
    SystemOfBodies bodies = createTestBodies( );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    const double observationTime = 5.0E8 + 1000.0;
    std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< ObservationModelSettings >( one_way_range, linkEnds ), bodies );

    std::vector< double > linkEndTimes;
    std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;
    const Eigen::VectorXd rangeObservation =
            observationModel->computeObservationsWithLinkEndData( observationTime, receiver, linkEndTimes, linkEndStates );

    BOOST_REQUIRE_EQUAL( linkEndTimes.size( ), 2 );
    const double expectedSolvedLightTime = rangeObservation( 0 ) / physical_constants::SPEED_OF_LIGHT;

    std::shared_ptr< ObservationDependentVariableSettings > lightTimeSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, receiver );
    std::shared_ptr< ObservationCollection<> > collection =
            simulateWithDependentVariables( one_way_range, bodies, linkEnds, { observationTime }, { lightTimeSettings } );

    const double computedLightTime = getScalarDependentVariable( collection, lightTimeSettings, observationTime );
    BOOST_CHECK_CLOSE_FRACTION( computedLightTime, expectedSolvedLightTime, 1.0e-13 );
    BOOST_CHECK_GT( std::fabs( ( linkEndTimes.at( 1 ) - linkEndTimes.at( 0 ) ) - computedLightTime ), 1.0e-10 );
}

//! Two-way uplink, downlink and combined propagation light times come from the per-leg solutions.
//! Combined light time excludes the retransmission delay.
BOOST_AUTO_TEST_CASE( testTwoWayLightTimeDependentVariableWithRetransmissionDelay )
{
    loadStandardSpiceKernels( );
    SystemOfBodies bodies = createTestBodies( );

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    twoWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    twoWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );

    const double observationTime = 5.0E8 + 1000.0;
    const double retransmissionDelay = 1.0e-3;
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings = getTwoWayRangeAncillarySettings( retransmissionDelay );

    std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< ObservationModelSettings >( n_way_range, twoWayLinkEnds ), bodies );

    std::vector< double > linkEndTimes;
    std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;
    const Eigen::VectorXd rangeObservation = observationModel->computeObservationsWithLinkEndData(
            observationTime, receiver, linkEndTimes, linkEndStates, ancillarySettings );

    BOOST_REQUIRE_EQUAL( linkEndTimes.size( ), 4 );
    const double expectedCombined = rangeObservation( 0 ) / physical_constants::SPEED_OF_LIGHT - retransmissionDelay;
    const double elapsedTransmitToReceive = linkEndTimes.at( 3 ) - linkEndTimes.at( 0 );
    const double observedRetransmissionDelay = linkEndTimes.at( 2 ) - linkEndTimes.at( 1 );

    BOOST_CHECK_SMALL( observedRetransmissionDelay - retransmissionDelay, 1.0e-6 );
    BOOST_CHECK_GT( std::fabs( elapsedTransmitToReceive - expectedCombined ), 0.5 * retransmissionDelay );

    std::shared_ptr< ObservationDependentVariableSettings > uplinkSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, retransmitter );
    std::shared_ptr< ObservationDependentVariableSettings > downlinkSettings =
            lightTimeBetweenLinkEndsDependentVariable( retransmitter, receiver );
    std::shared_ptr< ObservationDependentVariableSettings > combinedSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, receiver );
    std::shared_ptr< ObservationDependentVariableSettings > delaySettings = retransmissionDelaysDependentVariable( );
    std::shared_ptr< ObservationDependentVariableSettings > epochSettings = linkEndEpochsDependentVariable( );

    std::shared_ptr< ObservationCollection<> > collection =
            simulateWithDependentVariables( n_way_range,
                                            bodies,
                                            twoWayLinkEnds,
                                            { observationTime },
                                            { uplinkSettings, downlinkSettings, combinedSettings, delaySettings, epochSettings },
                                            ancillarySettings );

    const double uplinkLightTime = getScalarDependentVariable( collection, uplinkSettings, observationTime );
    const double downlinkLightTime = getScalarDependentVariable( collection, downlinkSettings, observationTime );
    const double combinedLightTime = getScalarDependentVariable( collection, combinedSettings, observationTime );

    BOOST_CHECK_CLOSE_FRACTION( combinedLightTime, expectedCombined, 1.0e-13 );
    BOOST_CHECK_CLOSE_FRACTION( combinedLightTime, uplinkLightTime + downlinkLightTime, 1.0e-13 );

    std::map< double, Eigen::VectorXd > delayHistory = collection->getDependentVariableHistory( delaySettings );
    BOOST_REQUIRE( delayHistory.count( observationTime ) == 1 );
    BOOST_REQUIRE_EQUAL( delayHistory.at( observationTime ).size( ), 1 );
    BOOST_CHECK_CLOSE_FRACTION(
            delayHistory.at( observationTime )( 0 ), retransmissionDelay, 10.0 * std::numeric_limits< double >::epsilon( ) );

    std::map< double, Eigen::VectorXd > epochHistory = collection->getDependentVariableHistory( epochSettings );
    BOOST_REQUIRE( epochHistory.count( observationTime ) == 1 );
    BOOST_REQUIRE_EQUAL( epochHistory.at( observationTime ).size( ), 4 );
    const double rangeAsTime = rangeObservation( 0 ) / physical_constants::SPEED_OF_LIGHT;
    BOOST_CHECK_CLOSE_FRACTION( rangeAsTime, combinedLightTime + retransmissionDelay, 1.0e-13 );
    BOOST_CHECK_CLOSE_FRACTION( rangeAsTime, elapsedTransmitToReceive, 1.0e-8 );
}

//! A three-leg path checks the generic role-to-leg mapping and verifies that every link-end
//! delay, including endpoint delays, is absent from the saved propagation-light-time total.
BOOST_AUTO_TEST_CASE( testThreeLegLightTimeExcludesAllLinkEndDelays )
{
    loadStandardSpiceKernels( );
    SystemOfBodies bodies = createTestBodies( );

    LinkEnds threeLegLinkEnds;
    threeLegLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    threeLegLinkEnds[ reflector1 ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    threeLegLinkEnds[ reflector2 ] = LinkEndId( std::make_pair( "Moon", "" ) );
    threeLegLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );

    const double observationTime = 5.0E8 + 1000.0;
    const std::vector< double > linkEndDelays{ 4.0e-4, 1.0e-3, 2.0e-3, 7.0e-4 };
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            getNWayRangeAncillarySettings( linkEndDelays );

    std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< ObservationModelSettings >( n_way_range, threeLegLinkEnds ), bodies );
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;
    const Eigen::VectorXd rangeObservation = observationModel->computeObservationsWithLinkEndData(
            observationTime, receiver, linkEndTimes, linkEndStates, ancillarySettings );

    BOOST_REQUIRE_EQUAL( linkEndTimes.size( ), 6 );
    const double totalLinkEndDelay = std::accumulate( linkEndDelays.begin( ), linkEndDelays.end( ), 0.0 );
    const double expectedPropagationLightTime =
            rangeObservation( 0 ) / physical_constants::SPEED_OF_LIGHT - totalLinkEndDelay;

    std::shared_ptr< ObservationDependentVariableSettings > firstLegSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, reflector1 );
    std::shared_ptr< ObservationDependentVariableSettings > middleLegSettings =
            lightTimeBetweenLinkEndsDependentVariable( reflector1, reflector2 );
    std::shared_ptr< ObservationDependentVariableSettings > lastLegSettings =
            lightTimeBetweenLinkEndsDependentVariable( reflector2, receiver );
    std::shared_ptr< ObservationDependentVariableSettings > firstTwoLegsSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, reflector2 );
    std::shared_ptr< ObservationDependentVariableSettings > fullPathSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, receiver );

    std::shared_ptr< ObservationCollection<> > collection = simulateWithDependentVariables(
            n_way_range,
            bodies,
            threeLegLinkEnds,
            { observationTime },
            { firstLegSettings, middleLegSettings, lastLegSettings, firstTwoLegsSettings, fullPathSettings },
            ancillarySettings );

    const double firstLeg = getScalarDependentVariable( collection, firstLegSettings, observationTime );
    const double middleLeg = getScalarDependentVariable( collection, middleLegSettings, observationTime );
    const double lastLeg = getScalarDependentVariable( collection, lastLegSettings, observationTime );
    const double firstTwoLegs = getScalarDependentVariable( collection, firstTwoLegsSettings, observationTime );
    const double fullPath = getScalarDependentVariable( collection, fullPathSettings, observationTime );

    BOOST_CHECK_CLOSE_FRACTION( fullPath, expectedPropagationLightTime, 1.0e-12 );
    BOOST_CHECK_CLOSE_FRACTION( firstTwoLegs, firstLeg + middleLeg, 1.0e-13 );
    BOOST_CHECK_CLOSE_FRACTION( fullPath, firstLeg + middleLeg + lastLeg, 1.0e-13 );
}

//! Integrated observables retain one solved leg cache per interval endpoint. Verify that the
//! dependent variable selects the requested cache instead of interpreting the doubled epoch array.
BOOST_AUTO_TEST_CASE( testIntegratedLightTimeSelectsRequestedIntervalEndpoint )
{
    loadStandardSpiceKernels( );
    SystemOfBodies bodies = createTestBodies( );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    const double observationTime = 5.0E8 + 1000.0;
    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings = getAveragedDopplerAncillarySettings( 60.0 );
    std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< ObservationModelSettings >( one_way_differenced_range, linkEnds ), bodies );
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;
    observationModel->computeObservationsWithLinkEndData(
            observationTime, receiver, linkEndTimes, linkEndStates, ancillarySettings );

    const auto legLightTimeCalculators = observationModel->getLegLightTimeCalculators( );
    const auto calculatorIt = legLightTimeCalculators.find( std::make_pair( transmitter, receiver ) );
    BOOST_REQUIRE( calculatorIt != legLightTimeCalculators.end( ) );
    BOOST_REQUIRE_EQUAL( calculatorIt->second.size( ), 2 );
    const double expectedStartLightTime = calculatorIt->second.at( 0 )->getCurrentTotalLightTime( );
    const double expectedEndLightTime = calculatorIt->second.at( 1 )->getCurrentTotalLightTime( );
    BOOST_CHECK_GT( std::fabs( expectedEndLightTime - expectedStartLightTime ), 1.0e-10 );

    std::shared_ptr< ObservationDependentVariableSettings > startSettings = lightTimeBetweenLinkEndsDependentVariable(
            transmitter, receiver, LinkEndId( "", "" ), LinkEndId( "", "" ), interval_start );
    std::shared_ptr< ObservationDependentVariableSettings > endSettings = lightTimeBetweenLinkEndsDependentVariable(
            transmitter, receiver, LinkEndId( "", "" ), LinkEndId( "", "" ), interval_end );
    std::shared_ptr< ObservationCollection<> > collection = simulateWithDependentVariables(
            one_way_differenced_range, bodies, linkEnds, { observationTime }, { startSettings, endSettings }, ancillarySettings );

    const double computedStartLightTime = getScalarDependentVariable( collection, startSettings, observationTime );
    const double computedEndLightTime = getScalarDependentVariable( collection, endSettings, observationTime );
    BOOST_CHECK_CLOSE_FRACTION( computedStartLightTime, expectedStartLightTime, 1.0e-13 );
    BOOST_CHECK_CLOSE_FRACTION( computedEndLightTime, expectedEndLightTime, 1.0e-13 );
}

//! An unspecified light-time request is a single, unambiguous full transmitter-to-receiver path.
BOOST_AUTO_TEST_CASE( testUnspecifiedLightTimeSelectsFullPath )
{
    loadStandardSpiceKernels( );
    SystemOfBodies bodies = createTestBodies( );

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    twoWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    twoWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );

    const double observationTime = 1.0E7 + 1000.0;
    std::shared_ptr< ObservationCollection<> > collection = simulateWithDependentVariables(
            n_way_range, bodies, twoWayLinkEnds, { observationTime }, { lightTimeBetweenLinkEndsDependentVariable( ) } );

    std::shared_ptr< SingleObservationSet< double, double > > singleSet =
            collection->getSingleLinkAndTypeObservationSets( n_way_range, LinkDefinition( twoWayLinkEnds ) ).at( 0 );
    BOOST_REQUIRE_EQUAL( singleSet->getDependentVariableBookkeeping( )->getDependentVariableSettings( ).size( ), 1 );

    const double defaultLightTime =
            getScalarDependentVariable( collection, lightTimeBetweenLinkEndsDependentVariable( ), observationTime );
    const double explicitFullPathLightTime = getScalarDependentVariable(
            collection, lightTimeBetweenLinkEndsDependentVariable( transmitter, receiver ), observationTime );
    BOOST_CHECK_CLOSE_FRACTION( defaultLightTime, explicitFullPathLightTime, 10.0 * std::numeric_limits< double >::epsilon( ) );
}

//! Exact ID matching lets a request for one station be applied across a batch that also contains
//! another station. A body-origin selector must not act as a wildcard for station endpoints.
BOOST_AUTO_TEST_CASE( testLightTimeSelectionUsesExactLinkEndIds )
{
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", "Station1" );
    linkEnds[ receiver ] = LinkEndId( "Earth", "Station1" );

    std::shared_ptr< ObservationDependentVariableSettings > station1Settings = lightTimeBetweenLinkEndsDependentVariable(
            transmitter, receiver, LinkEndId( "Earth", "Station1" ), LinkEndId( "Earth", "Station1" ) );
    std::shared_ptr< ObservationDependentVariableSettings > station2Settings = lightTimeBetweenLinkEndsDependentVariable(
            transmitter, receiver, LinkEndId( "Earth", "Station2" ), LinkEndId( "Earth", "Station2" ) );
    std::shared_ptr< ObservationDependentVariableSettings > bodyOriginSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, receiver, LinkEndId( "Earth" ), LinkEndId( "Earth" ) );

    BOOST_REQUIRE_EQUAL( createAllCompatibleDependentVariableSettings( one_way_range, linkEnds, station1Settings ).size( ), 1 );
    BOOST_CHECK( createAllCompatibleDependentVariableSettings( one_way_range, linkEnds, station2Settings ).empty( ) );
    BOOST_CHECK( createAllCompatibleDependentVariableSettings( one_way_range, linkEnds, bodyOriginSettings ).empty( ) );

    std::shared_ptr< ObservationDependentVariableSettings > startOnlySettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter );
    const std::vector< std::shared_ptr< ObservationDependentVariableSettings > > startOnlyCompleted =
            createAllCompatibleDependentVariableSettings( one_way_range, linkEnds, startOnlySettings );
    BOOST_REQUIRE_EQUAL( startOnlyCompleted.size( ), 1 );
    const std::shared_ptr< InterlinkObservationDependentVariableSettings > completedInterlinkSettings =
            std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( startOnlyCompleted.at( 0 ) );
    BOOST_REQUIRE( completedInterlinkSettings != nullptr );
    BOOST_CHECK_EQUAL( completedInterlinkSettings->originatingLinkEndType_, transmitter );
    BOOST_CHECK_EQUAL( completedInterlinkSettings->linkEndType_, receiver );
}

//! Reversed start/end types on a two-way observable are rejected explicitly.
BOOST_AUTO_TEST_CASE( testReversedLightTimeSelectionThrows )
{
    loadStandardSpiceKernels( );
    SystemOfBodies bodies = createTestBodies( );

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    twoWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    twoWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );

    std::shared_ptr< TabulatedObservationSimulationSettings<> > simulationSettings =
            std::make_shared< TabulatedObservationSimulationSettings<> >( n_way_range, twoWayLinkEnds, std::vector< double >{ 1.0E7 } );
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput{ simulationSettings };

    BOOST_CHECK_THROW( addDependentVariablesToObservationSimulationSettings(
                               measurementSimulationInput, { lightTimeBetweenLinkEndsDependentVariable( receiver, transmitter ) }, bodies ),
                       std::runtime_error );
}

//! Solved one-way light time with a Shapiro correction matches the observation model's solved
//! light time, not the uncorrected solution.
BOOST_AUTO_TEST_CASE( testOneWayLightTimeIncludesConfiguredCorrections )
{
    loadStandardSpiceKernels( );
    SystemOfBodies bodies = createTestBodies( );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    const double observationTime = 5.0E8 + 1000.0;
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > correctionSettings{
        std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Sun" } )
    };

    std::shared_ptr< ObservationModel< 1, double, double > > uncorrectedModel =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< ObservationModelSettings >( one_way_range, linkEnds ), bodies );
    std::shared_ptr< ObservationModel< 1, double, double > > correctedModel =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< ObservationModelSettings >( one_way_range, linkEnds, correctionSettings ), bodies );

    std::vector< double > uncorrectedTimes;
    std::vector< Eigen::Matrix< double, 6, 1 > > uncorrectedStates;
    const Eigen::VectorXd uncorrectedRange =
            uncorrectedModel->computeObservationsWithLinkEndData( observationTime, receiver, uncorrectedTimes, uncorrectedStates );

    std::vector< double > correctedTimes;
    std::vector< Eigen::Matrix< double, 6, 1 > > correctedStates;
    const Eigen::VectorXd correctedRange =
            correctedModel->computeObservationsWithLinkEndData( observationTime, receiver, correctedTimes, correctedStates );

    const double uncorrectedLightTime = uncorrectedRange( 0 ) / physical_constants::SPEED_OF_LIGHT;
    const double correctedLightTime = correctedRange( 0 ) / physical_constants::SPEED_OF_LIGHT;
    BOOST_CHECK_GT( std::fabs( correctedLightTime - uncorrectedLightTime ), 1.0e-12 );

    std::shared_ptr< ObservationDependentVariableSettings > lightTimeSettings =
            lightTimeBetweenLinkEndsDependentVariable( transmitter, receiver );
    std::shared_ptr< ObservationCollection<> > collection = simulateWithDependentVariables(
            one_way_range, bodies, linkEnds, { observationTime }, { lightTimeSettings }, nullptr, correctionSettings );

    const double computedLightTime = getScalarDependentVariable( collection, lightTimeSettings, observationTime );
    BOOST_CHECK_CLOSE_FRACTION( computedLightTime, correctedLightTime, 1.0e-13 );
    BOOST_CHECK_GT( std::fabs( computedLightTime - uncorrectedLightTime ), 1.0e-12 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
