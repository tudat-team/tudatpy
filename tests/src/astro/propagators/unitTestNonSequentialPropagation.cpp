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

    // Check consistency of state history sizes
    BOOST_CHECK_EQUAL( nonsequentialStateHistory.size( ), forwardStateHistory.size( ) + backwardStateHistory.size( ) - 1 );
    BOOST_CHECK_EQUAL( nonsequentialDependentVariablesHistory.size( ), forwardDependentVariablesHistory.size( ) + backwardDependentVariablesHistory.size( ) - 1 );

    // Check propagated state and dependent variables at initial time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonsequentialStateHistory.begin( )->second, backwardStateHistory.begin( )->second, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonsequentialDependentVariablesHistory.begin( )->second, backwardDependentVariablesHistory.begin( )->second,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check propagated state and dependent variables at final time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonsequentialStateHistory.rbegin( )->second, forwardStateHistory.rbegin( )->second, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonsequentialDependentVariablesHistory.rbegin( )->second, forwardDependentVariablesHistory.rbegin( )->second,
                                       std::numeric_limits< double >::epsilon( ) );

}


BOOST_AUTO_TEST_CASE( testNonSequentialMultiArcDynamics )
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
    std::vector< double > arcStartTimes = { 0.0, 8.0 * 3600.0, 16.0 * 3600.0 };
    std::vector< double > arcEndTimes = { 8.0 * 3600.0, 16.0 * 3600.0, 86400.0 };
    std::vector< double > midArcTimes;
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        midArcTimes.push_back( ( arcStartTimes.at( i ) + arcEndTimes.at( i ) ) / 2.0 );
    }

    unsigned int nbArcs = arcStartTimes.size( );

    std::string globalFrameOrientation = "ECLIPJ2000";
    std::string globalFrameOrigin = "Jupiter";

    // Create bodies.
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEpoch - 86400.0, finalEpoch + 86400.0 );
    bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );

    bodySettings.at( "Io" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Europa" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Ganymede" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Callisto" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

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
    double timeStep = 2700.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > forwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > backwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, - timeStep, CoefficientSets::rungeKuttaFehlberg78, - timeStep, - timeStep, 1.0e3, 1.0e3 );

    // Define arc-wise initial states
    std::vector< Eigen::VectorXd > midArcStatesMoons;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        midArcStatesMoons.push_back( propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcTimes.at( i ) ) );
    }

    // Define dependent variables to save
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", "Jupiter" ) );
    for ( unsigned int k = 0 ; k < bodiesToPropagate.size( ) ; k++ )
    {
        if ( bodiesToPropagate.at( k ) != "Io" )
        {
            dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", bodiesToPropagate.at( k ) ) );
        }
    }


    // Define propagator settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > forwardPropagationSettingsList, backwardPropagationSettingsList, nonsequentialPropagationSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        forwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), arcEndTimes.at( i ), cowell, dependentVariables ) );
        backwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), arcStartTimes.at( i ), cowell, dependentVariables ) );

        std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings =
                std::make_shared< NonSequentialPropagationTerminationSettings >(
                        std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ),
                        std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ) );
        std::shared_ptr< TranslationalStatePropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), forwardIntegratorSettings,
                terminationSettings, cowell, dependentVariables );
        nonsequentialPropagationSettingsList.push_back( nonsequentialPropagatorSettings );
    }
    std::shared_ptr< MultiArcPropagatorSettings< > > forwardPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( forwardPropagationSettingsList );
    std::shared_ptr< MultiArcPropagatorSettings< > > backwardPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( backwardPropagationSettingsList );

    std::shared_ptr< MultiArcPropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( nonsequentialPropagationSettingsList );

    // Propagate dynamics (forward leg)
    bool setIntegratedResult = true;
    MultiArcDynamicsSimulator< > forwardDynamicsSimulator = MultiArcDynamicsSimulator< >(
            bodies, forwardIntegratorSettings, forwardPropagatorSettings, midArcTimes, true, false, setIntegratedResult );

    // Propagate dynamics (backward leg)
    MultiArcDynamicsSimulator< > backwardDynamicsSimulator = MultiArcDynamicsSimulator< >(
            bodies, backwardIntegratorSettings, backwardPropagatorSettings, midArcTimes, true, false, setIntegratedResult );

    // Save propagations outputs
    std::vector< std::map< double, Eigen::VectorXd > > forwardStateHistory = forwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::vector< std::map< double, Eigen::VectorXd > > forwardDependentVariablesHistory = forwardDynamicsSimulator.getDependentVariableHistory( );
    std::vector< std::map< double, Eigen::VectorXd > > backwardStateHistory = backwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::vector< std::map< double, Eigen::VectorXd > > backwardDependentVariablesHistory = backwardDynamicsSimulator.getDependentVariableHistory( );

    // Propagate dynamics (non-sequentially)
    MultiArcDynamicsSimulator< > nonSequentialDynamicsSimulator = MultiArcDynamicsSimulator< >( bodies, nonsequentialPropagatorSettings, true );

    // Save non-sequential propagations outputs
    std::vector< std::map< double, Eigen::VectorXd > > nonsequentialStateHistory = nonSequentialDynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::vector< std::map< double, Eigen::VectorXd > > nonsequentialDependentVariablesHistory = nonSequentialDynamicsSimulator.getDependentVariableHistory( );

    for ( unsigned int i = 0 ; i < nonsequentialStateHistory.size( ) ; i++ )
    {
        BOOST_CHECK_EQUAL( nonsequentialStateHistory.at( i ).size( ), forwardStateHistory.at( i ).size( ) + backwardStateHistory.at( i ).size( ) - 1 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonsequentialStateHistory.at( i ).begin( )->second, backwardStateHistory.at( i ).begin( )->second,
                                           std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonsequentialStateHistory.at( i ).rbegin( )->second, forwardStateHistory.at( i ).rbegin( )->second,
                                           std::numeric_limits< double >::epsilon( ) );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}