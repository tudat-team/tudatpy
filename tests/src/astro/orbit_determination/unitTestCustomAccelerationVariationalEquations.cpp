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

#include <string>
#include <thread>

#include <limits>


#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_custom_acceleration_estimation )

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::ground_stations;
using namespace tudat::observation_models;


Eigen::Vector3d customAccelerationFunction(
    const SystemOfBodies& bodies,
    const double time )
{
    double scaleDistance = 1.0E6;
    Eigen::Vector6d relativeState = bodies.at( "Vehicle" )->getState( ) - bodies.at( "Earth" )->getState( );
    Eigen::Vector3d acceleration = Eigen::Vector3d::Ones( ) * std::exp( -relativeState.segment( 0, 3 ).norm( ) / scaleDistance );
    return acceleration;
}

Eigen::MatrixXd getCustomAccelerationPartialFunction(
    const SystemOfBodies& bodies,
    const double time,
    const Eigen::Vector3d& acceleration )
{
    double scaleDistance = 1.0E6;
    Eigen::Vector6d relativeState = bodies.at( "Vehicle" )->getState( ) - bodies.at( "Earth" )->getState( );
    Eigen::Matrix< double, 3, 6 > partial = Eigen::Matrix< double, 3, 6 >::Zero( );
    double positionNorm = relativeState.segment( 0, 3 ).norm( );
    for( unsigned int i = 0; i < 3; i ++ )
    {
        partial.block( i, 0, 1, 3 ) = -relativeState.segment( 0, 3 ).transpose( ) / ( positionNorm * scaleDistance ) *
            std::exp( -positionNorm / scaleDistance );
    }

    return partial;
}

BOOST_AUTO_TEST_CASE( test_CustomAccelerationEstimation )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );
    double finalEphemerisTime = initialEphemerisTime + 1.0 * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< CustomAccelerationSettings >(
        std::bind( &customAccelerationFunction, bodies, std::placeholders::_1 ) ) );
//    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
//                                                     basic_astrodynamics::point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
//    dependentVariables.push_back(
//                std::make_shared< SingleDependentVariableSaveSettings >(
//                    aerodynamic_force_coefficients_dependent_variable, "Vehicle", "Earth" ) );
//    dependentVariables.push_back(
//                std::make_shared< SingleDependentVariableSaveSettings >(
//                    radiation_pressure_coefficient_dependent_variable, "Vehicle", "Sun" ) );


    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
        rungeKuttaFixedStepSettings( 90.0, CoefficientSets::rungeKuttaFehlberg78 );

    for( int testCase = 0; testCase < 3; testCase ++ )
    {

        std::shared_ptr<TranslationalStatePropagatorSettings<double> > propagatorSettings =
            std::make_shared<TranslationalStatePropagatorSettings<double> >(
                centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime,
                integratorSettings,
                std::make_shared<PropagationTimeTerminationSettings>( finalEphemerisTime ),
                cowell, dependentVariables );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Define list of parameters to estimate.
        std::vector<std::shared_ptr<EstimatableParameterSettings> > parameterNames;
        parameterNames = getInitialStateParameterSettings<double>( propagatorSettings, bodies );
        if( testCase == 1 )
        {
            parameterNames.at( 0 )->customPartialSettings_ =
                std::make_shared<NumericalAccelerationPartialSettings>( Eigen::Vector6d::Ones( ) * 10.0,
                                                                        "Vehicle", "Vehicle", custom_acceleration );
        }
        else if( testCase == 2 )
        {
            parameterNames.at( 0 )->customPartialSettings_ =
                std::make_shared<AnalyticalAccelerationPartialSettings>(
                    std::bind( &getCustomAccelerationPartialFunction, bodies, std::placeholders::_1, std::placeholders::_2 ),
                    "Vehicle", "Vehicle", custom_acceleration );
        }

        // Create parameters
        std::shared_ptr<estimatable_parameters::EstimatableParameterSet<double> > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////          INITIALIZE ORBIT DETERMINATION OBJECT     ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        SingleArcVariationalEquationsSolver<> variationalEquationsSolver = SingleArcVariationalEquationsSolver<>(
            bodies, propagatorSettings, parametersToEstimate );

        std::map<double, Eigen::MatrixXd>
            stateTransitionMatrixHistory = variationalEquationsSolver.getStateTransitionMatrixSolution( );

        double positionPerturbation = 10.0;
        double velocityPerturbation = 1.0E-2;

        double statePerturbation;
        Eigen::VectorXd unperturbedInitialState = propagatorSettings->getInitialStates( );
        Eigen::VectorXd perturbedInitialState;
        std::map<double, Eigen::VectorXd> upperturbedStateHistory, downperturbedStateHistory;
        std::map<double, Eigen::MatrixXd> numericalStateTransitionMatrixHistory;

        for ( unsigned int i = 0; i < 6; i++ )
        {
            if ( i < 3 )
            {
                statePerturbation = positionPerturbation;
            }
            else
            {
                statePerturbation = velocityPerturbation;
            }

            perturbedInitialState = unperturbedInitialState;
            perturbedInitialState( i ) += statePerturbation;

            propagatorSettings->resetInitialStates( perturbedInitialState );
            SingleArcDynamicsSimulator<> upperturbedDynamicsSimulator( bodies, propagatorSettings );
            upperturbedStateHistory = upperturbedDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            perturbedInitialState = unperturbedInitialState;
            perturbedInitialState( i ) -= statePerturbation;

            propagatorSettings->resetInitialStates( perturbedInitialState );
            SingleArcDynamicsSimulator<> downperturbedDynamicsSimulator( bodies, propagatorSettings );
            downperturbedStateHistory = downperturbedDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            for ( auto it: upperturbedStateHistory )
            {
                if ( downperturbedStateHistory.count( it.first ) == 0 )
                {
                    throw std::runtime_error( "Error in custom acceleration partial test, state size is inconsistent" );
                }
                if ( stateTransitionMatrixHistory.count( it.first ) == 0 )
                {
                    throw std::runtime_error(
                        "Error in custom acceleration partial test, variational equations size is inconsistent" );
                }
                if ( i == 0 )
                {
                    numericalStateTransitionMatrixHistory[ it.first ] = Eigen::Matrix6d::Zero( );
                }

                numericalStateTransitionMatrixHistory[ it.first ].block( 0, i, 6, 1 ) =
                    ( upperturbedStateHistory.at( it.first ) - downperturbedStateHistory.at( it.first )) /
                    ( 2.0 * statePerturbation );

            }
        }


        double maximumError = 0.0;
        for ( auto it: numericalStateTransitionMatrixHistory )
        {
            auto matrix1 = it.second;
            auto matrix2 = stateTransitionMatrixHistory.at( it.first );

            for ( int j = 0; j < 6; j++ )
            {
                int rowBlock = ( j < 3 ) ? 0 : 1;
                for ( int k = 0; k < 6; k++ )
                {
                    int columnBlock = ( k < 3 ) ? 0 : 1;
                    double entryDifference = matrix1( j, k ) - matrix2( j, k );
                    double blockNorm = matrix1.block( rowBlock * 3, columnBlock * 3, 3, 3 ).norm( );
//                    BOOST_CHECK_SMALL( std::fabs( entryDifference ), blockNorm * 1.0E-6 );

                    if ( std::fabs( entryDifference ) / blockNorm > maximumError )
                    {
                        maximumError = std::fabs( entryDifference ) / blockNorm;
                    }
                }
            }
        }
        if( testCase == 0 )
        {
            BOOST_CHECK_SMALL( 1.0 / maximumError, 50.0 );
        }
        else if( testCase == 1 )
        {
            BOOST_CHECK_SMALL( maximumError, 1.0E-6 );
        }
        else if( testCase == 2 )
        {
            BOOST_CHECK_SMALL( maximumError, 1.0E-6 );
        }
        std::cout << maximumError << std::endl;
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
