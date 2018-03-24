/* git    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"

namespace tudat
{

namespace unit_tests
{

//Using declarations.
using namespace tudat;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( test_variational_equation_calculation )


template< typename TimeType = double , typename StateScalarType  = double >
        std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
executeMultiArcEarthMoonSimulation(
        const std::vector< std::string > centralBodies,
        const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference = Eigen::Matrix< StateScalarType, 12, 1 >::Zero( ),
        const int propagationType = 0,
        const bool patchArcsTogether = 0,
        const Eigen::Vector3d parameterPerturbation = Eigen::Vector3d::Zero( ),
        const bool propagateVariationalEquations = 1,
        const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > forcedArcInitialStates =
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( ),
        const double arcDuration = 5.0E5,
        const double arcOverlap  = 5.0E3 )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = initialEphemerisTime + 0.5E7;
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings[ "Moon" ]->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings[ "Earth" ]->ephemerisSettings->resetMakeMultiArcEphemeris( true );


    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Earth" ] = accelerationsOfEarth;

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    accelerationsOfMoon[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMoon[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Moon" ] = accelerationsOfMoon;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Moon" );
    bodiesToIntegrate.push_back( "Earth" );

    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define propagator settings.
    std::map< std::string, std::string > centralBodyMap;

    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, centralBodyMap );

    // Create integrator settings
    boost::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            boost::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime ), 1800.0 );


    // Define arc times.
    std::vector< double > arcStartTimes, arcEndTimes;
    TimeType integrationStartTime = initialEphemerisTime;
    TimeType integrationEndTime = finalEphemerisTime;

    TimeType currentStartTime = integrationStartTime;
    TimeType currentEndTime = integrationStartTime + arcDuration;

    do
    {
        arcStartTimes.push_back( currentStartTime );
        arcEndTimes.push_back( currentEndTime );

        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < ( integrationEndTime ) );


    // Set (perturbed) initial state.
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > systemInitialStates;
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        if( forcedArcInitialStates.size( ) == 0 )
        {
            systemInitialStates.push_back( getInitialStatesOfBodies< TimeType, StateScalarType >(
                                               bodiesToIntegrate, centralBodies, bodyMap, arcStartTimes.at( i ) ) );
        }
        else
        {
            systemInitialStates.push_back( forcedArcInitialStates.at( i ) );
        }
        systemInitialStates[ i ] += initialStateDifference;
    }

    // Create propagator settings
    TranslationalPropagatorType propagatorType;
    if( propagationType == 0 )
    {
        propagatorType = cowell;
    }
    else if( propagationType == 1 )
    {
        propagatorType = encke;
    }
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsList;

    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        propagatorSettingsList.push_back( boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                                          ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialStates.at( i ),
                                            arcEndTimes.at( i ), propagatorType ) );
    }
    boost::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcPropagatorSettings;

    if( patchArcsTogether && forcedArcInitialStates.size( ) == 0 )
    {
        multiArcPropagatorSettings =
                boost::make_shared< MultiArcPropagatorSettings< StateScalarType > >(
                    propagatorSettingsList, patchArcsTogether );
    }
    else
    {
        multiArcPropagatorSettings =
                boost::make_shared< MultiArcPropagatorSettings< StateScalarType > >(
                    propagatorSettingsList );
    }

    // Define parameters.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames.push_back(
                    boost::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        "Moon", arcStartTimes, centralBodies[ 0 ] ) );
        parameterNames.push_back(
                    boost::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        "Earth", arcStartTimes, centralBodies[ 1 ] ) );
        parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );
        parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter ) );
        parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );

    }

    // Create parameters
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Perturb parameters.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );

    parameterVector.block( parameterVector.rows( ) - 3, 0, 3, 1 ) += parameterPerturbation;
    parametersToEstimate->resetParameterValues( parameterVector );

    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;
    {
        // Create dynamics simulator
        MultiArcVariationalEquationsSolver< StateScalarType, TimeType > variationalEquations =
                MultiArcVariationalEquationsSolver< StateScalarType, TimeType >(
                    bodyMap, integratorSettings, multiArcPropagatorSettings, parametersToEstimate, arcStartTimes );

        // Propagate requested equations.
        if( propagateVariationalEquations )
        {
            variationalEquations.integrateVariationalAndDynamicalEquations( multiArcPropagatorSettings->getInitialStateList( ), 1 );
        }
        else
        {
            variationalEquations.integrateDynamicalEquationsOfMotionOnly( multiArcPropagatorSettings->getInitialStateList( ) );
        }


        // Retrieve test data
        for( unsigned int arc = 0; arc < arcEndTimes.size( ); arc++ )
        {
            double testEpoch = arcEndTimes.at( arc ) - 2.0E4;
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates =
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
            testStates.block( 0, 0, 6, 1 ) = bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( testEpoch );

            if( centralBodyMap[ "Moon" ] == "Earth" )
            {
                testStates.block( 0, 0, 6, 1 ) -= bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( testEpoch );
            }

            testStates.block( 6, 0, 6, 1 ) = bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( testEpoch );

            if( propagateVariationalEquations )
            {
                results.first.push_back( variationalEquations.getStateTransitionMatrixInterface( )->
                                         getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
                results.second.push_back( multiArcPropagatorSettings->getInitialStateList( ).at( arc ) );
            }
            else
            {
                results.second.push_back( testStates );
            }
        }
    }
    return results;
}

BOOST_AUTO_TEST_CASE( testEarthMoonMultiArcVariationalEquationCalculation )
{
    std::pair< std::vector< Eigen::MatrixXd >, std::vector< Eigen::VectorXd > > currentOutput;

    std::vector< std::vector< std::string > > centralBodiesSet;
    std::vector< std::string > centralBodies;

    // Define central bodt settings
    centralBodies.resize( 2 );

    centralBodies[ 0 ] = "SSB";
    centralBodies[ 1 ] = "SSB";
    centralBodiesSet.push_back( centralBodies );

    centralBodies[ 0 ] = "Earth";
    centralBodies[ 1 ] = "Sun";
    centralBodiesSet.push_back( centralBodies );


    // Define variables for numerical differentiation
    Eigen::Matrix< double, 12, 1>  perturbedState;
    Eigen::Vector3d perturbedParameter;

    Eigen::Matrix< double, 12, 1> statePerturbation;
    Eigen::Vector3d parameterPerturbation;


    for( unsigned int i = 0; i < 2; i++ )
    {
        // Define parameter perturbation
        parameterPerturbation  = ( Eigen::Vector3d( ) << 1.0E11, 1.0E11, 1.0E15 ).finished( );

        // Define state perturbation
        if( i == 0 )
        {
            statePerturbation = ( Eigen::Matrix< double, 12, 1>( ) <<
                                  100000.0, 100000.0, 100000.0, 1.0, 0.1, 1.0,
                                  100000.0, 100000.0, 100000.0, 1.0, 0.1, 1.0 ).finished( );
        }
        else if( i == 1 )
        {
            statePerturbation = ( Eigen::Matrix< double, 12, 1>( ) <<
                                  100000.0, 100000.0, 100000.0, 1.0, 0.1, 1.0,
                                  10000000.0, 1000000.0, 100000000.0, 100.0, 100.0, 1000.0 ).finished( );
        }

        // Only perform Encke propagation is origin is not SSB
        unsigned int maximumPropagatorType = 1;
        if( i == 1 )
        {
            maximumPropagatorType = 1;
        }

        for( unsigned int patchArcs = 0; patchArcs < 2; patchArcs++ )
        {
            // Test for all requested propagator types.
            for( unsigned int k = 0; k < maximumPropagatorType; k++ )
            {
                std::cout << "Test: " << i << " " << patchArcs << " " << k << std::endl;
                // Compute state transition and sensitivity matrices
                currentOutput = executeMultiArcEarthMoonSimulation < double, double >(
                            centralBodiesSet[ i ], Eigen::Matrix< double, 12, 1 >::Zero( ), k, patchArcs );
                std::vector< Eigen::MatrixXd > stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first;
                std::vector< Eigen::VectorXd > nominalArcStartStates;
                if( patchArcs )
                {
                    nominalArcStartStates = currentOutput.second;
                }

                std::vector< Eigen::MatrixXd > manualPartial;
                manualPartial.resize( stateTransitionAndSensitivityMatrixAtEpoch.size( ) );
                for( unsigned int arc = 0; arc < manualPartial.size( ); arc++ )
                {
                    manualPartial[ arc ] = Eigen::MatrixXd::Zero( 12, 15 );
                }
                // Numerically compute state transition matrix
                for( unsigned int j = 0; j < 12; j++ )
                {
                    std::vector< Eigen::VectorXd > upPerturbedState, downPerturbedState;
                    perturbedState.setZero( );
                    perturbedState( j ) += statePerturbation( j );
                    upPerturbedState = executeMultiArcEarthMoonSimulation< double, double >(
                                centralBodiesSet[ i ], perturbedState, k, patchArcs, Eigen::Vector3d::Zero( ), 0, nominalArcStartStates ).second;

                    perturbedState.setZero( );
                    perturbedState( j ) -= statePerturbation( j );
                    downPerturbedState = executeMultiArcEarthMoonSimulation< double, double >(
                                centralBodiesSet[ i ], perturbedState, k, patchArcs, Eigen::Vector3d::Zero( ), 0, nominalArcStartStates ).second;

                    for( unsigned int arc = 0; arc < upPerturbedState.size( ); arc++ )
                    {
                        manualPartial[ arc ].block( 0, j, 12, 1 ) =
                                ( upPerturbedState[ arc ] - downPerturbedState[ arc ] ) / ( 2.0 * statePerturbation( j ) );
                    }
                }

                // Numerically compute sensitivity matrix
                for( unsigned int j = 0; j < 3; j ++ )
                {
                    std::vector< Eigen::VectorXd > upPerturbedState, downPerturbedState;
                    perturbedState.setZero( );
                    Eigen::Vector3d upPerturbedParameter, downPerturbedParameter;
                    perturbedParameter.setZero( );
                    perturbedParameter( j ) += parameterPerturbation( j );
                    upPerturbedState = executeMultiArcEarthMoonSimulation< double, double >(
                                centralBodiesSet[ i ], perturbedState, k, patchArcs, perturbedParameter, 0, nominalArcStartStates ).second;

                    perturbedParameter.setZero( );
                    perturbedParameter( j ) -= parameterPerturbation( j );
                    downPerturbedState = executeMultiArcEarthMoonSimulation< double, double >(
                                centralBodiesSet[ i ], perturbedState, k, patchArcs, perturbedParameter, 0, nominalArcStartStates ).second;

                    for( unsigned int arc = 0; arc < upPerturbedState.size( ); arc++ )
                    {
                        manualPartial[ arc ].block( 0, j + 12, 12, 1 ) =
                                ( upPerturbedState[ arc ] - downPerturbedState[ arc ] ) / ( 2.0 * parameterPerturbation( j ) );
                    }
                }

                std::cout << manualPartial.at( 1 ) << std::endl << std::endl
                          << stateTransitionAndSensitivityMatrixAtEpoch.at( 1 ) << std::endl;

                // Check results
                for( unsigned int arc = 0; arc < manualPartial.size( ); arc++ )
                {
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 0, 0, 12, 15 ),
                                manualPartial.at( arc ).block( 0, 0, 12, 15 ), 5.0E-4 );
                }

            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}

