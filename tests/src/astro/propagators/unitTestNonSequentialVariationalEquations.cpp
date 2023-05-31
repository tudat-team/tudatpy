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
#include "tudat/simulation/estimation.h"
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
using namespace estimatable_parameters;

//! Test suite for astro functions.
BOOST_AUTO_TEST_SUITE( test_non_sequential_variational_equations )

BOOST_AUTO_TEST_CASE( testNonSequentialSingleArcVariationalEquations )
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
    std::shared_ptr< TranslationalStatePropagatorSettings< > > forwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, finalEpoch, cowell, dependentVariables );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > backwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, initialEpoch, cowell, dependentVariables );

    // Define parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back(
                std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                        bodiesToPropagate.at( i ), midArcStatesMoons.segment( i * 6, 6 ), centralBodies.at( i ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                bodiesToPropagate.at( i ), gravitational_parameter ) );

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
    }

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > forwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, forwardPropagatorSettings );
    printEstimatableParameterEntries( forwardPropagationParametersToEstimate );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > backwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, backwardPropagatorSettings );

    // Propagate variational equations (forward leg)
    bool setIntegratedResult = false;
    bool integrateConcurrently = true;
    SingleArcVariationalEquationsSolver< > forwardVariationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
            bodies, forwardIntegratorSettings, forwardPropagatorSettings, forwardPropagationParametersToEstimate, integrateConcurrently,
            forwardIntegratorSettings, false, true, setIntegratedResult );

    // Propagate variational equations (backward leg)
    SingleArcVariationalEquationsSolver< > backwardVariationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
            bodies, backwardIntegratorSettings, backwardPropagatorSettings, backwardPropagationParametersToEstimate, integrateConcurrently,
            backwardIntegratorSettings, false, true, setIntegratedResult );

    // Save propagations outputs
    std::map< double, Eigen::MatrixXd > forwardStmHistory = forwardVariationalEquationsSolver.getStateTransitionMatrixSolution( );
    std::map< double, Eigen::MatrixXd > forwardSemHistory = forwardVariationalEquationsSolver.getSensitivityMatrixSolution( );
    std::map< double, Eigen::MatrixXd > backwardStmHistory = backwardVariationalEquationsSolver.getStateTransitionMatrixSolution( );
    std::map< double, Eigen::MatrixXd > backwardSemHistory = backwardVariationalEquationsSolver.getSensitivityMatrixSolution( );

    //! Create settings for non-sequential propagation
    std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
            std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ), std::make_shared< PropagationTimeTerminationSettings >( initialEpoch ) );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, terminationSettings, cowell, dependentVariables );


    // Propagate variational equations (non-sequentially)
    SingleArcVariationalEquationsSolver< > nonsequentialVariationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
            bodies, forwardIntegratorSettings, nonsequentialPropagatorSettings, forwardPropagationParametersToEstimate, integrateConcurrently,
            forwardIntegratorSettings, false, true, setIntegratedResult );

    std::map< double, Eigen::MatrixXd > nonsequentialStmHistory = nonsequentialVariationalEquationsSolver.getStateTransitionMatrixSolution( );
    std::map< double, Eigen::MatrixXd > nonsequentialSemHistory = nonsequentialVariationalEquationsSolver.getSensitivityMatrixSolution( );

    // Check consistency maps size
    BOOST_CHECK_EQUAL( nonsequentialStmHistory.size( ), forwardStmHistory.size( ) + backwardStmHistory.size( ) - 1 );
    BOOST_CHECK_EQUAL( nonsequentialSemHistory.size( ), forwardSemHistory.size( ) + backwardStmHistory.size( ) - 1 );

    // Check state transition and sensitivity matrices at initial time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( backwardStmHistory.begin( )->second, nonsequentialStmHistory.begin( )->second, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( backwardSemHistory.begin( )->second, nonsequentialSemHistory.begin( )->second, std::numeric_limits< double >::epsilon( ) );

    // Check state transition and sensitivity matrices at final time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( forwardStmHistory.rbegin( )->second, nonsequentialStmHistory.rbegin( )->second, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( forwardSemHistory.rbegin( )->second, nonsequentialSemHistory.rbegin( )->second, std::numeric_limits< double >::epsilon( ) );


}

BOOST_AUTO_TEST_SUITE_END( )

}

}