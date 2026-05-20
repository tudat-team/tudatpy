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

#include <limits>
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

BOOST_AUTO_TEST_SUITE( test_light_time_correction_components )

//! Helper: simulate a scalar range observable (one-way or n-way) at a list of times,
//! with the given light-time corrections, and (optionally) record the per-correction contributions
//! as a dependent variable. Returns the observation history and the dependent-variable history.
struct SimulationOutputs
{
    std::map< double, Eigen::VectorXd > observations;
    std::map< double, Eigen::VectorXd > dependentVariables;
};

SimulationOutputs simulateRangeObservableWithDependentVariables(
        const ObservableType observableType,
        const SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const std::vector< double >& observationTimes,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& correctionSettings,
        std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesList,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
        const LinkEndType referenceLinkEnd = receiver )
{
    std::shared_ptr< ObservationModelSettings > observationSettings;
    if( observableType == n_way_differenced_range )
    {
        observationSettings =
                std::make_shared< NWayDifferencedRangeObservationModelSettings >( linkEnds, correctionSettings, nullptr );
    }
    else
    {
        observationSettings = std::make_shared< ObservationModelSettings >( observableType, linkEnds, correctionSettings );
    }
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
            createObservationSimulators( std::vector< std::shared_ptr< ObservationModelSettings > >{ observationSettings }, bodies );

    std::shared_ptr< TabulatedObservationSimulationSettings<> > simulationSettings =
            std::make_shared< TabulatedObservationSimulationSettings<> >(
                    observableType, linkEnds, observationTimes, referenceLinkEnd, std::vector< std::shared_ptr< ObservationViabilitySettings > >{ },
                    nullptr, ancillarySettings );
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput{ simulationSettings };

    if( dependentVariablesList.empty( ) == false )
    {
        addDependentVariablesToObservationSimulationSettings( measurementSimulationInput, dependentVariablesList, bodies );
    }

    std::shared_ptr< ObservationCollection<> > collection =
            simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );
    std::shared_ptr< SingleObservationSet< double, double > > singleSet =
            collection->getSingleLinkAndTypeObservationSets( observableType, LinkDefinition( linkEnds ) ).at( 0 );

    SimulationOutputs outputs;
    outputs.observations = singleSet->getObservationsHistory( );
    if( dependentVariablesList.empty( ) == false )
    {
        outputs.dependentVariables = singleSet->getDependentVariableHistory( );
    }
    return outputs;
}

SimulationOutputs simulateRangeObservable(
        const ObservableType observableType,
        const SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const std::vector< double >& observationTimes,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& correctionSettings,
        const std::shared_ptr< ObservationDependentVariableSettings >& dependentVariableSettings = nullptr )
{
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesList;
    if( dependentVariableSettings != nullptr )
    {
        dependentVariablesList.push_back( dependentVariableSettings );
    }

    return simulateRangeObservableWithDependentVariables(
            observableType, bodies, linkEnds, observationTimes, correctionSettings, dependentVariablesList );
}

//! Verifies that clearing dependent-variable settings also clears deferred settings whose size
//! has not yet been resolved from an observation model.
BOOST_AUTO_TEST_CASE( testClearSettingsClearsDeferredLightTimeComponents )
{
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    std::shared_ptr< ObservationDependentVariableBookkeeping > bookkeeping =
            std::make_shared< ObservationDependentVariableBookkeeping >( one_way_range, LinkDefinition( linkEnds ) );
    bookkeeping->addDependentVariables( { lightTimeCorrectionComponentsDependentVariable( transmitter, receiver ) } );
    BOOST_REQUIRE_EQUAL( bookkeeping->getDeferredSettings( ).size( ), 1 );

    bookkeeping->clearSettings( );
    BOOST_CHECK( bookkeeping->getDeferredSettings( ).empty( ) );
    BOOST_CHECK_EQUAL( bookkeeping->getDependentVariableSettings( ).size( ), 0 );
    BOOST_CHECK_EQUAL( bookkeeping->getTotalDependentVariableSize( ), 0 );
}

//! Verifies that an unspecified light-time correction component request expands to the direct
//! observable leg, not to the reversed leg.
BOOST_AUTO_TEST_CASE( testUnspecifiedLightTimeComponentLegUsesDirectLeg )
{
    loadStandardSpiceKernels( );

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
    createGroundStation(
            bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    std::shared_ptr< LightTimeCorrectionSettings > sunCorrection =
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Sun" } );

    SimulationOutputs outputs = simulateRangeObservable(
            one_way_range,
            bodies,
            linkEnds,
            { 1.0E7 + 1000.0 },
            { sunCorrection },
            lightTimeCorrectionComponentsDependentVariable( ) );

    BOOST_REQUIRE_EQUAL( outputs.dependentVariables.size( ), 1 );
    BOOST_REQUIRE_EQUAL( outputs.dependentVariables.begin( )->second.size( ), 1 );
    BOOST_CHECK_GT( outputs.dependentVariables.begin( )->second( 0 ), 0.0 );
}

//! Verifies that the `light_time_correction_components` dependent variable on a one-way range
//! observable reproduces the individual light-time correction contributions. Two relativistic
//! corrections are used (perturbers: Sun and Jupiter), so the vector has two entries.
//!
//! Four simulations are compared:
//!  1. No corrections.
//!  2. Sun correction only.
//!  3. Jupiter correction only.
//!  4. Both corrections — also saves the per-correction dependent variable.
//!
//! Checks (per observation epoch):
//!   - The dependent-variable vector has length 2 and is non-NaN/finite.
//!   - Each component is positive (first-order relativistic Shapiro delays on a light-like path
//!     are positive semi-definite for masses at a distance from the path).
//!   - `c × component_i ≈ range(only_i) − range(none)` — each saved component reproduces the
//!     effect of enabling just that correction.
//!   - `c × (component_0 + component_1) ≈ range(both) − range(none)` — the sum reproduces the
//!     effect of enabling both.
//!
//! The tolerance accounts for the second-order residual introduced by light-time iteration
//! re-converging to slightly different transmission times when the correction set changes. The
//! residual scales as (correction × velocity_of_link_end), which for Earth-station → Moon-orbiter
//! with solar Shapiro delay is well below a millimetre at typical geometries.
BOOST_AUTO_TEST_CASE( testPerCorrectionComponentsOneWayRange )
{
    loadStandardSpiceKernels( );

    // Set up the environment. MoonOrbiter is placed in a simple Keplerian orbit around the Moon.
    std::vector< std::string > bodyNames{ "Earth", "Moon", "Sun", "Jupiter" };
    const double initialEphemerisTime = 1.0E7;
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth" );

    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings =
            keplerEphemerisSettings( keplerElements, 0.0, getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    createGroundStation(
            bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );

    // One-way range: Earth-Station1 (transmitter) → MoonOrbiter (receiver).
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    std::vector< double > observationTimes;
    for( int i = 0; i < 10; i++ )
    {
        observationTimes.push_back( initialEphemerisTime + 1000.0 + i * 100.0 );
    }

    // Two relativistic corrections, one per perturber.
    std::shared_ptr< LightTimeCorrectionSettings > sunCorrection =
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Sun" } );
    std::shared_ptr< LightTimeCorrectionSettings > jupiterCorrection =
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Jupiter" } );

    // 1. No corrections.
    SimulationOutputs noCorrections = simulateRangeObservable( one_way_range, bodies, linkEnds, observationTimes, { } );
    // 2. Sun correction only.
    SimulationOutputs sunOnly = simulateRangeObservable( one_way_range, bodies, linkEnds, observationTimes, { sunCorrection } );
    // 3. Jupiter correction only.
    SimulationOutputs jupiterOnly = simulateRangeObservable( one_way_range, bodies, linkEnds, observationTimes, { jupiterCorrection } );
    // 4. Both corrections + dependent variable that saves per-correction components.
    SimulationOutputs both = simulateRangeObservable( one_way_range,
                                                      bodies,
                                                      linkEnds,
                                                      observationTimes,
                                                      { sunCorrection, jupiterCorrection },
                                                      lightTimeCorrectionComponentsDependentVariable( transmitter, receiver ) );
    // 5. Sun correction only + dependent variable.
    SimulationOutputs sunOnlyWithDependent = simulateRangeObservable(
            one_way_range, bodies, linkEnds, observationTimes, { sunCorrection }, lightTimeCorrectionComponentsDependentVariable( transmitter, receiver ) );
    // 6. Jupiter correction only + dependent variable.
    SimulationOutputs jupiterOnlyWithDependent = simulateRangeObservable( one_way_range,
                                                                          bodies,
                                                                          linkEnds,
                                                                          observationTimes,
                                                                          { jupiterCorrection },
                                                                          lightTimeCorrectionComponentsDependentVariable( transmitter, receiver ) );

    BOOST_REQUIRE_EQUAL( noCorrections.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( sunOnly.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( jupiterOnly.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( both.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( both.dependentVariables.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( sunOnlyWithDependent.dependentVariables.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( jupiterOnlyWithDependent.dependentVariables.size( ), observationTimes.size( ) );

    const double c = physical_constants::SPEED_OF_LIGHT;

    // Tolerance on the residual between (range_with_correction - range_without) and c × component.
    // The residual is dominated by the second-order change in converged transmission time when
    // corrections are enabled; ~(correction × link-end velocity). 1 mm is a comfortable bound
    // for Earth-station → Moon-orbiter with relativistic corrections only.
    const double rangeTolerance = 1.0E-3;

    for( auto it: both.dependentVariables )
    {
        const double time = it.first;
        const Eigen::VectorXd& components = it.second;

        // The dependent-variable vector holds one entry per registered correction.
        BOOST_REQUIRE_EQUAL( components.size( ), 2 );

        // Values are finite (not NaN).
        BOOST_CHECK( components( 0 ) == components( 0 ) );
        BOOST_CHECK( components( 1 ) == components( 1 ) );

        // First-order-relativistic Shapiro delay is positive for a mass off the line of sight.
        BOOST_CHECK_GT( components( 0 ), 0.0 );
        BOOST_CHECK_GT( components( 1 ), 0.0 );

        // Sun's Shapiro delay on an Earth–Moon leg dominates over Jupiter's.
        BOOST_CHECK_GT( components( 0 ), components( 1 ) );

        // Per-correction equivalence.
        const double rangeNone = noCorrections.observations.at( time )( 0 );
        const double rangeSun = sunOnly.observations.at( time )( 0 );
        const double rangeJupiter = jupiterOnly.observations.at( time )( 0 );
        const double rangeBoth = both.observations.at( time )( 0 );

        BOOST_CHECK_SMALL( std::fabs( c * components( 0 ) - ( rangeSun - rangeNone ) ), rangeTolerance );
        BOOST_CHECK_SMALL( std::fabs( c * components( 1 ) - ( rangeJupiter - rangeNone ) ), rangeTolerance );

        // Sum equivalence: sum of saved components reproduces the effect of enabling both corrections.
        BOOST_CHECK_SMALL(
                std::fabs( c * ( components( 0 ) + components( 1 ) ) - ( rangeBoth - rangeNone ) ), rangeTolerance );

        // Single-correction saves produce 1-entry vectors equal to their isolated range effect.
        const Eigen::VectorXd& sunComponentOnly = sunOnlyWithDependent.dependentVariables.at( time );
        const Eigen::VectorXd& jupiterComponentOnly = jupiterOnlyWithDependent.dependentVariables.at( time );
        BOOST_REQUIRE_EQUAL( sunComponentOnly.size( ), 1 );
        BOOST_REQUIRE_EQUAL( jupiterComponentOnly.size( ), 1 );
        BOOST_CHECK_SMALL( std::fabs( c * sunComponentOnly( 0 ) - ( rangeSun - rangeNone ) ), rangeTolerance );
        BOOST_CHECK_SMALL( std::fabs( c * jupiterComponentOnly( 0 ) - ( rangeJupiter - rangeNone ) ), rangeTolerance );
        BOOST_CHECK_CLOSE_FRACTION( sunComponentOnly( 0 ), components( 0 ), std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION(
                jupiterComponentOnly( 0 ), components( 1 ), std::numeric_limits< double >::epsilon( ) );
    }
}

//! Verifies per-leg `light_time_correction_components` behaviour for a two-leg n-way range
//! observable (transmitter -> retransmitter -> transmitter) with a single correction model.
BOOST_AUTO_TEST_CASE( testPerCorrectionComponentsTwoWayNWayRange )
{
    loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames{ "Earth", "Moon", "Sun" };
    const double initialEphemerisTime = 1.0E7;
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth" );

    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings =
            keplerEphemerisSettings( keplerElements, 0.0, getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    createGroundStation(
            bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    twoWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    twoWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );

    std::vector< double > observationTimes{ initialEphemerisTime + 1000.0,
                                            initialEphemerisTime + 1100.0,
                                            initialEphemerisTime + 1200.0 };

    std::shared_ptr< LightTimeCorrectionSettings > sunCorrection =
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Sun" } );

    SimulationOutputs noCorrections =
            simulateRangeObservable( n_way_range, bodies, twoWayLinkEnds, observationTimes, { } );
    SimulationOutputs sunOnly =
            simulateRangeObservable( n_way_range, bodies, twoWayLinkEnds, observationTimes, { sunCorrection } );

    SimulationOutputs uplinkOnly = simulateRangeObservable(
            n_way_range,
            bodies,
            twoWayLinkEnds,
            observationTimes,
            { sunCorrection },
            lightTimeCorrectionComponentsDependentVariable(
                    transmitter,
                    retransmitter,
                    LinkEndId( std::make_pair( "Earth", "Station1" ) ),
                    LinkEndId( std::make_pair( "MoonOrbiter", "" ) ) ) );

    SimulationOutputs downlinkOnly = simulateRangeObservable(
            n_way_range,
            bodies,
            twoWayLinkEnds,
            observationTimes,
            { sunCorrection },
            lightTimeCorrectionComponentsDependentVariable(
                    retransmitter,
                    receiver,
                    LinkEndId( std::make_pair( "MoonOrbiter", "" ) ),
                    LinkEndId( std::make_pair( "Earth", "Station1" ) ) ) );

    BOOST_REQUIRE_EQUAL( noCorrections.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( sunOnly.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( uplinkOnly.dependentVariables.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( downlinkOnly.dependentVariables.size( ), observationTimes.size( ) );

    const double c = physical_constants::SPEED_OF_LIGHT;
    const double rangeTolerance = 1.0E-3;

    for( const auto& entry : uplinkOnly.dependentVariables )
    {
        const double time = entry.first;
        const Eigen::VectorXd& upComponent = entry.second;
        const Eigen::VectorXd& downComponent = downlinkOnly.dependentVariables.at( time );

        BOOST_REQUIRE_EQUAL( upComponent.size( ), 1 );
        BOOST_REQUIRE_EQUAL( downComponent.size( ), 1 );

        const double rangeNone = noCorrections.observations.at( time )( 0 );
        const double rangeSun = sunOnly.observations.at( time )( 0 );
        BOOST_CHECK_SMALL( std::fabs( c * ( upComponent( 0 ) + downComponent( 0 ) ) - ( rangeSun - rangeNone ) ), rangeTolerance );
    }
}

//! Verifies that each leg of an n-way differenced range exposes the light-time correction
//! components from both its arc-start and arc-end range calculators.
BOOST_AUTO_TEST_CASE( testPerCorrectionComponentsNWayDifferencedRange )
{
    loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames{ "Earth", "Moon", "Sun" };
    const double initialEphemerisTime = 1.0E7;
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth" );

    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings =
            keplerEphemerisSettings( keplerElements, 0.0, getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    createGroundStation(
            bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );

    LinkEndId earthStation( std::make_pair( "Earth", "Station1" ) );
    LinkEndId moonOrbiter( std::make_pair( "MoonOrbiter", "" ) );

    LinkEnds twoWayLinkEnds;
    twoWayLinkEnds[ transmitter ] = earthStation;
    twoWayLinkEnds[ retransmitter ] = moonOrbiter;
    twoWayLinkEnds[ receiver ] = earthStation;

    const double integrationTime = 60.0;
    std::vector< double > observationTimes{ initialEphemerisTime + 1000.0,
                                            initialEphemerisTime + 1100.0,
                                            initialEphemerisTime + 1200.0 };

    std::vector< double > rangeObservationTimes;
    for( double observationTime: observationTimes )
    {
        rangeObservationTimes.push_back( observationTime - integrationTime / 2.0 );
        rangeObservationTimes.push_back( observationTime + integrationTime / 2.0 );
    }

    std::shared_ptr< LightTimeCorrectionSettings > sunCorrection =
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Sun" } );

    SimulationOutputs rangeNoCorrections =
            simulateRangeObservable( n_way_range, bodies, twoWayLinkEnds, rangeObservationTimes, { } );
    SimulationOutputs rangeSunOnly =
            simulateRangeObservable( n_way_range, bodies, twoWayLinkEnds, rangeObservationTimes, { sunCorrection } );

    SimulationOutputs differencedNoCorrections = simulateRangeObservableWithDependentVariables(
            n_way_differenced_range,
            bodies,
            twoWayLinkEnds,
            observationTimes,
            { },
            std::vector< std::shared_ptr< ObservationDependentVariableSettings > >{ },
            getAveragedDopplerAncillarySettings( integrationTime ) );

    SimulationOutputs differencedSunOnly = simulateRangeObservableWithDependentVariables(
            n_way_differenced_range,
            bodies,
            twoWayLinkEnds,
            observationTimes,
            { sunCorrection },
            std::vector< std::shared_ptr< ObservationDependentVariableSettings > >{
                    lightTimeCorrectionComponentsDependentVariable( transmitter, retransmitter, earthStation, moonOrbiter ),
                    lightTimeCorrectionComponentsDependentVariable( retransmitter, receiver, moonOrbiter, earthStation ) },
            getAveragedDopplerAncillarySettings( integrationTime ) );

    BOOST_REQUIRE_EQUAL( rangeNoCorrections.observations.size( ), rangeObservationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( rangeSunOnly.observations.size( ), rangeObservationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( differencedNoCorrections.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( differencedSunOnly.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( differencedSunOnly.dependentVariables.size( ), observationTimes.size( ) );

    const double c = physical_constants::SPEED_OF_LIGHT;
    const double rangeTolerance = 1.0E-3;
    const double differencedRangeTolerance = 1.0E-4;

    for( const auto& entry: differencedSunOnly.dependentVariables )
    {
        const double time = entry.first;
        const Eigen::VectorXd& components = entry.second;

        BOOST_REQUIRE_EQUAL( components.size( ), 4 );
        for( int i = 0; i < components.size( ); i++ )
        {
            BOOST_CHECK( components( i ) == components( i ) );
            BOOST_CHECK_GT( components( i ), 0.0 );
        }

        const double startTime = time - integrationTime / 2.0;
        const double endTime = time + integrationTime / 2.0;
        const double startRangeCorrection =
                rangeSunOnly.observations.at( startTime )( 0 ) - rangeNoCorrections.observations.at( startTime )( 0 );
        const double endRangeCorrection =
                rangeSunOnly.observations.at( endTime )( 0 ) - rangeNoCorrections.observations.at( endTime )( 0 );

        const double startDependentCorrection = c * ( components( 0 ) + components( 2 ) );
        const double endDependentCorrection = c * ( components( 1 ) + components( 3 ) );
        BOOST_CHECK_SMALL( std::fabs( startDependentCorrection - startRangeCorrection ), rangeTolerance );
        BOOST_CHECK_SMALL( std::fabs( endDependentCorrection - endRangeCorrection ), rangeTolerance );

        const double differencedCorrection =
                differencedSunOnly.observations.at( time )( 0 ) - differencedNoCorrections.observations.at( time )( 0 );
        const double dependentDifferencedCorrection =
                ( endDependentCorrection - startDependentCorrection ) / integrationTime;
        BOOST_CHECK_SMALL( std::fabs( dependentDifferencedCorrection - differencedCorrection ), differencedRangeTolerance );
    }
}

//! Verifies that the `correction_type_filter` argument selects a subset of the registered
//! corrections by `LightTimeCorrectionType`. Two first-order-relativistic corrections with
//! different perturbers are registered, then filtered by type.
BOOST_AUTO_TEST_CASE( testCorrectionTypeFilter )
{
    loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames{ "Earth", "Moon", "Sun", "Jupiter" };
    const double initialEphemerisTime = 1.0E7;
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth" );

    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings =
            keplerEphemerisSettings( keplerElements, 0.0, getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    createGroundStation(
            bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    std::vector< double > observationTimes{ initialEphemerisTime + 1000.0, initialEphemerisTime + 1100.0 };

    std::shared_ptr< LightTimeCorrectionSettings > sunCorrection =
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Sun" } );
    std::shared_ptr< LightTimeCorrectionSettings > jupiterCorrection =
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Jupiter" } );

    // Full vector (no filter): 2 entries.
    SimulationOutputs full = simulateRangeObservable( one_way_range,
                                                      bodies,
                                                      linkEnds,
                                                      observationTimes,
                                                      { sunCorrection, jupiterCorrection },
                                                      lightTimeCorrectionComponentsDependentVariable( transmitter, receiver ) );

    // Filtered (first-order relativistic only): all selected entries must match the unfiltered
    // vector for this setup, because all registered corrections are first-order relativistic.
    SimulationOutputs filtered = simulateRangeObservable(
            one_way_range,
            bodies,
            linkEnds,
            observationTimes,
            { sunCorrection, jupiterCorrection },
            lightTimeCorrectionComponentsDependentVariable(
                    transmitter,
                    receiver,
                    LinkEndId( "", "" ),
                    LinkEndId( "", "" ),
                    std::vector< LightTimeCorrectionType >{ first_order_relativistic } ) );

    BOOST_REQUIRE_EQUAL( full.dependentVariables.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( filtered.dependentVariables.size( ), observationTimes.size( ) );

    for( auto it: filtered.dependentVariables )
    {
        const double time = it.first;
        const Eigen::VectorXd& fullComponents = full.dependentVariables.at( time );
        BOOST_REQUIRE_EQUAL( it.second.size( ), fullComponents.size( ) );
        for( int i = 0; i < it.second.size( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( it.second( i ), fullComponents( i ), std::numeric_limits< double >::epsilon( ) );
        }
    }
}

//! Verifies that angular position exposes its light-time leg through the generic observation-model
//! accessor and can therefore save `light_time_correction_components`.
BOOST_AUTO_TEST_CASE( testAngularPositionCorrectionComponents )
{
    loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames{ "Earth", "Moon", "Sun" };
    const double initialEphemerisTime = 1.0E7;
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth" );

    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings =
            keplerEphemerisSettings( keplerElements, 0.0, getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    createGroundStation(
            bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", "Station1" ) );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > correctionSettings{
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >{ "Sun" } )
    };
    std::shared_ptr< ObservationModelSettings > observationSettings = std::make_shared< ObservationModelSettings >(
            angular_position, linkEnds, correctionSettings );
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
            createObservationSimulators( std::vector< std::shared_ptr< ObservationModelSettings > >{ observationSettings }, bodies );

    std::shared_ptr< TabulatedObservationSimulationSettings<> > simulationSettings =
            std::make_shared< TabulatedObservationSimulationSettings<> >(
                    angular_position, linkEnds, std::vector< double >{ initialEphemerisTime + 1000.0 }, receiver );
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput{ simulationSettings };

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesList{
            lightTimeCorrectionComponentsDependentVariable( transmitter, receiver )
    };
    addDependentVariablesToObservationSimulationSettings( measurementSimulationInput, dependentVariablesList, bodies );

    std::shared_ptr< ObservationCollection<> > collection =
            simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );
    std::shared_ptr< SingleObservationSet< double, double > > singleSet =
            collection->getSingleLinkAndTypeObservationSets( angular_position, LinkDefinition( linkEnds ) ).at( 0 );

    std::map< double, Eigen::VectorXd > dependentVariableHistory = singleSet->getDependentVariableHistory( );
    BOOST_REQUIRE_EQUAL( dependentVariableHistory.size( ), 1 );
    BOOST_REQUIRE_EQUAL( dependentVariableHistory.begin( )->second.size( ), 1 );
    BOOST_CHECK( dependentVariableHistory.begin( )->second( 0 ) == dependentVariableHistory.begin( )->second( 0 ) );
    BOOST_CHECK_GT( dependentVariableHistory.begin( )->second( 0 ), 0.0 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
