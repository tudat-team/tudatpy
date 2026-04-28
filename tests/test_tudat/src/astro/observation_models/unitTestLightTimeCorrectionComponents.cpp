/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
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

//! Helper: simulate a one-way range observable for Earth-Station1 → MoonOrbiter at a list of times,
//! with the given light-time corrections, and (optionally) record the per-correction contributions
//! as a dependent variable. Returns the observation history and the dependent-variable history.
struct SimulationOutputs
{
    std::map< double, Eigen::VectorXd > observations;
    std::map< double, Eigen::VectorXd > dependentVariables;
};

SimulationOutputs simulateOneWayRange(
        const SystemOfBodies& bodies,
        const LinkEnds& linkEnds,
        const std::vector< double >& observationTimes,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& correctionSettings,
        const std::shared_ptr< ObservationDependentVariableSettings >& dependentVariableSettings = nullptr )
{
    std::shared_ptr< ObservationModelSettings > observationSettings =
            std::make_shared< ObservationModelSettings >( one_way_range, linkEnds, correctionSettings );
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
            createObservationSimulators( std::vector< std::shared_ptr< ObservationModelSettings > >{ observationSettings }, bodies );

    std::shared_ptr< TabulatedObservationSimulationSettings<> > simulationSettings =
            std::make_shared< TabulatedObservationSimulationSettings<> >( one_way_range, linkEnds, observationTimes, receiver );
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput{ simulationSettings };

    if( dependentVariableSettings != nullptr )
    {
        std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesList{ dependentVariableSettings };
        addDependentVariablesToObservationSimulationSettings( measurementSimulationInput, dependentVariablesList, bodies );
    }

    std::shared_ptr< ObservationCollection<> > collection =
            simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );
    std::shared_ptr< SingleObservationSet< double, double > > singleSet =
            collection->getSingleLinkAndTypeObservationSets( one_way_range, LinkDefinition( linkEnds ) ).at( 0 );

    SimulationOutputs outputs;
    outputs.observations = singleSet->getObservationsHistory( );
    if( dependentVariableSettings != nullptr )
    {
        outputs.dependentVariables = singleSet->getDependentVariableHistory( );
    }
    return outputs;
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
    SimulationOutputs noCorrections = simulateOneWayRange( bodies, linkEnds, observationTimes, { } );
    // 2. Sun correction only.
    SimulationOutputs sunOnly = simulateOneWayRange( bodies, linkEnds, observationTimes, { sunCorrection } );
    // 3. Jupiter correction only.
    SimulationOutputs jupiterOnly = simulateOneWayRange( bodies, linkEnds, observationTimes, { jupiterCorrection } );
    // 4. Both corrections + dependent variable that saves per-correction components.
    SimulationOutputs both = simulateOneWayRange(
            bodies,
            linkEnds,
            observationTimes,
            { sunCorrection, jupiterCorrection },
            lightTimeCorrectionComponentsDependentVariable( transmitter, receiver ) );

    BOOST_REQUIRE_EQUAL( noCorrections.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( sunOnly.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( jupiterOnly.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( both.observations.size( ), observationTimes.size( ) );
    BOOST_REQUIRE_EQUAL( both.dependentVariables.size( ), observationTimes.size( ) );

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
    }
}

//! Verifies that the `correction_type_filter` argument selects a subset of the registered
//! corrections by `LightTimeCorrectionType`. Two corrections of the same type (first-order
//! relativistic, but with different perturbers) are registered, then filtered — the filter
//! matches by type, so the result has one entry corresponding to the first-registered correction
//! of that type.
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
    SimulationOutputs full = simulateOneWayRange(
            bodies,
            linkEnds,
            observationTimes,
            { sunCorrection, jupiterCorrection },
            lightTimeCorrectionComponentsDependentVariable( transmitter, receiver ) );

    // Filtered (first-order relativistic only): the filter matches by type; with two corrections
    // of the same type, the factory picks the first match, yielding a 1-entry vector equal to the
    // first component of the full vector.
    SimulationOutputs filtered = simulateOneWayRange(
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
        BOOST_REQUIRE_EQUAL( it.second.size( ), 1 );
        BOOST_CHECK_CLOSE_FRACTION(
                it.second( 0 ), full.dependentVariables.at( time )( 0 ), std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
