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

#include "tudat/simulation/simulation.h"
#include <boost/test/unit_test.hpp>
#include "tudat/basics/testMacros.h"

namespace tudat
{
namespace unit_tests
{

//! Using declarations.
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace ephemerides;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;

//! Test suite for astro functions.
BOOST_AUTO_TEST_SUITE( test_non_sequential_propagation )

BOOST_AUTO_TEST_CASE( testNonSequentialSingleArcDynamics )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Callisto" );

    // Specify initial time
    double initialEpoch = 0.0;
    double finalEpoch = 86400.0;
    double midArcEpoch = ( initialEpoch + finalEpoch ) / 2.0;

    std::string globalFrameOrientation = "ECLIPJ2000";
    std::string globalFrameOrigin = "Jupiter";

    // Create bodies.
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEpoch - 86400.0, finalEpoch + 86400.0 );
    bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    //! Define bodies to propagate.
    std::vector< std::string > bodiesToPropagate = { "Io", "Europa", "Ganymede", "Callisto" };
    std::vector< std::string > centralBodies = { "Jupiter", "Jupiter", "Jupiter", "Jupiter" };

    // Set accelerations.
    SelectedAccelerationMap accelerationSettings;
    for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 8, 0, 2, 2  )  );
        for ( unsigned int j = 0 ; j < bodiesToPropagate.size( ) ; j++ )
        {
            if ( i != j )
            {
                accelerationsOfSatellite[ bodiesToPropagate[ j ] ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 8, 0 ) );
            }
        }
        accelerationSettings[ bodiesToPropagate[ i ] ] = accelerationsOfSatellite;
    }

    basic_astrodynamics::AccelerationMap accelerationsMap = createAccelerationModelsMap(
            bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    // Define integrator settings
    double initialTimeIntegration = midArcEpoch;
    double timeStep = 3600.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > forwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    initialTimeIntegration, timeStep, CoefficientSets::rungeKuttaFehlberg78,
                    timeStep, timeStep, 1.0e3, 1.0e3 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > backwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    initialTimeIntegration, - timeStep, CoefficientSets::rungeKuttaFehlberg78,
                    - timeStep, - timeStep, 1.0e3, 1.0e3 );

    // Define initial states
    Eigen::VectorXd midArcStatesMoons = propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcEpoch );

    // Define dependent variables to save
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesSettings;
    dependentVariablesSettings.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", "Jupiter" ) );
    for ( unsigned int k = 0 ; k < bodiesToPropagate.size( ) ; k++ )
    {
        if ( bodiesToPropagate.at( k ) != "Io" )
        {
            dependentVariablesSettings.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", bodiesToPropagate.at( k ) ) );
        }
    }

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< > > forwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, finalEpoch, cowell, dependentVariablesSettings );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > backwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, initialEpoch, cowell, dependentVariablesSettings );

    // Propagate dynamics (forward leg)
    bool setIntegratedResult = false;
    SingleArcDynamicsSimulator< > forwardDynamicsSimulator = SingleArcDynamicsSimulator< >(
            bodies, forwardIntegratorSettings, forwardPropagatorSettings, true, false, setIntegratedResult, true );

    // Propagate dynamics (backward leg)
    SingleArcDynamicsSimulator< > backwardDynamicsSimulator = SingleArcDynamicsSimulator< >(
            bodies, backwardIntegratorSettings, backwardPropagatorSettings, true, false, setIntegratedResult, true );

    // Save propagations outputs
    std::map< double, Eigen::VectorXd > forwardStateHistory = forwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > forwardDependentVariablesHistory = forwardDynamicsSimulator.getDependentVariableHistory( );
    std::map< double, Eigen::VectorXd > backwardStateHistory = backwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > backwardDependentVariablesHistory = backwardDynamicsSimulator.getDependentVariableHistory( );

    std::map< double, Eigen::VectorXd > combinedStateHistory = forwardStateHistory;
    combinedStateHistory.insert( backwardStateHistory.begin( ), backwardStateHistory.end( ) );
    std::map< double, Eigen::VectorXd > combinedDependentVariablesHistory = forwardDependentVariablesHistory;
    combinedDependentVariablesHistory.insert( backwardDependentVariablesHistory.begin( ), backwardDependentVariablesHistory.end( ) );

    //! Create settings for non-sequential propagation
    std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings =
            std::make_shared< NonSequentialPropagationTerminationSettings >(
                    std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ),
                    std::make_shared< PropagationTimeTerminationSettings >( initialEpoch ) );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, terminationSettings, cowell, dependentVariablesSettings );

    // Propagate dynamics (non-sequentially)
    SingleArcDynamicsSimulator< > nonsequentialDynamicsSimulator = SingleArcDynamicsSimulator< >(
            bodies, forwardIntegratorSettings, nonsequentialPropagatorSettings, true, false, true, true );

    std::map< double, Eigen::VectorXd > nonsequentialStateHistory = nonsequentialDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > nonsequentialDependentVariablesHistory = nonsequentialDynamicsSimulator.getDependentVariableHistory( );

    for ( auto itr : combinedStateHistory )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( itr.second, nonsequentialStateHistory.at( itr.first ), std::numeric_limits< double >::epsilon( ) );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}