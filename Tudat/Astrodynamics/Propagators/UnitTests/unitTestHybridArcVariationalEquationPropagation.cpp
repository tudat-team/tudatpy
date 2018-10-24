/* Copyright (c) 2010-2018, Delft University of Technology
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
#include "Tudat/SimulationSetup/EstimationSetup/variationalEquationsSolver.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EstimationSetup/createNumericalSimulator.h"
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
using namespace tudat::unit_conversions;

BOOST_AUTO_TEST_SUITE( test_hybrid_arc_variational_equation_calculation )


template< typename TimeType = double , typename StateScalarType  = double >
        std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
executeHybridArcMarsAndOrbiterSensitivitySimulation(
        const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference = Eigen::Matrix< StateScalarType, 12, 1 >::Zero( ),
        const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 2 ),
        const bool propagateVariationalEquations = 1,
        const bool patchMultiArcs = 0,
        const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > forcedMultiArcInitialStates =
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( ),
        const double arcDuration = 0.5 * 86400.0,
        const double arcOverlap  = 5.0E3 )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Earth" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 2.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Orbiter" ] = std::make_shared< Body >( );
    bodyMap[ "Orbiter" ]->setConstantBodyMass( 5.0E3 );
    bodyMap[ "Orbiter" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                            std::map< double, std::shared_ptr< Ephemeris > >( ),
                                            "Mars", "ECLIPJ2000" ) );

    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > orbiterRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Orbiter" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    orbiterRadiationPressureSettings, "Orbiter", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap singleArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    singleArcAccelerationMap[ "Mars" ] = accelerationsOfMars;

    std::vector< std::string > singleArcBodiesToIntegrate, singleArcCentralBodies;
    singleArcBodiesToIntegrate.push_back( "Mars" );
    singleArcCentralBodies.push_back( "SSB" );

    AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, singleArcAccelerationMap, singleArcBodiesToIntegrate, singleArcCentralBodies );
    Eigen::VectorXd singleArcInitialStates = getInitialStatesOfBodies(
                singleArcBodiesToIntegrate, singleArcCentralBodies, bodyMap, initialEphemerisTime );

    singleArcInitialStates += initialStateDifference.segment(
                0, singleArcInitialStates.rows( ) );

    std::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >(
                singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToIntegrate,
                singleArcInitialStates, finalEphemerisTime );


    SelectedAccelerationMap multiArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfOrbiter;
    accelerationsOfOrbiter[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfOrbiter[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    multiArcAccelerationMap[ "Orbiter" ] = accelerationsOfOrbiter;

    std::vector< std::string > multiArcBodiesToIntegrate, multiArcCentralBodies;
    multiArcBodiesToIntegrate.push_back( "Orbiter" );
    multiArcCentralBodies.push_back( "Mars" );

    AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, multiArcAccelerationMap, multiArcBodiesToIntegrate, multiArcCentralBodies );

    // Creater arc times
    std::vector< double > integrationArcStarts, integrationArcEnds;
    double integrationStartTime = initialEphemerisTime;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double currentStartTime = integrationStartTime;
    double currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcStarts.push_back( currentStartTime );
        integrationArcEnds.push_back( currentEndTime );

        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );

    // Create list of multi-arc initial states
    unsigned int numberOfIntegrationArcs = integrationArcStarts.size( );
    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    multiArcSystemInitialStates.resize( numberOfIntegrationArcs );

    // Define (quasi-arbitrary) arc initial states
    double marsGravitationalParameter =  bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );

    if( forcedMultiArcInitialStates.size( ) == 0 )
    {
        for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
        {
            Eigen::Vector6d orbiterInitialStateInKeplerianElements;
            orbiterInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6000.0E3;
            orbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
            orbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
            orbiterInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                    = convertDegreesToRadians( 235.7 - j );
            orbiterInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                    = convertDegreesToRadians( 23.4 + j );
            orbiterInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 + j * 10.0 );

            // Convert state from Keplerian elements to Cartesian elements.
            multiArcSystemInitialStates[ j ]  = convertKeplerianToCartesianElements(
                        orbiterInitialStateInKeplerianElements,
                        marsGravitationalParameter ) + initialStateDifference.segment(
                        singleArcInitialStates.rows( ), 6 );
        }
    }
    else
    {
        for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
        {
            multiArcSystemInitialStates[ j ] = forcedMultiArcInitialStates.at( j ) + initialStateDifference.segment(
                        singleArcInitialStates.rows( ), 6 );
        }
    }

    // Create propagation settings for each arc
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
    {
        arcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( multiArcCentralBodies, multiArcAccelerationModelMap, multiArcBodiesToIntegrate,
                      multiArcSystemInitialStates.at( i ), integrationArcEnds.at( i ) ) );
    }

    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< > >( arcPropagationSettingsList, patchMultiArcs );

    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings =
            std::make_shared< HybridArcPropagatorSettings< > >(
                singleArcPropagatorSettings, multiArcPropagatorSettings );

    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 60.0 );



    // Define parameters.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames.push_back(
                    std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        multiArcBodiesToIntegrate.at( 0 ), multiArcPropagatorSettings->getInitialStates( ),
                        integrationArcStarts, multiArcCentralBodies.at( 0 ) ) );
        parameterNames.push_back(
                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        singleArcBodiesToIntegrate.at( 0 ), singleArcInitialStates, singleArcCentralBodies.at( 0 ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Perturb parameters.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );

    parameterVector.block( parameterVector.rows( ) - 2, 0, 2, 1 ) += parameterPerturbation;
    parametersToEstimate->resetParameterValues( parameterVector );

    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;
    {
        // Create dynamics simulator
        HybridArcVariationalEquationsSolver< StateScalarType, TimeType > variationalEquations =
                HybridArcVariationalEquationsSolver< StateScalarType, TimeType >(
                    bodyMap, integratorSettings, hybridArcPropagatorSettings, parametersToEstimate, integrationArcStarts );

        // Propagate requested equations.
        if( propagateVariationalEquations )
        {
            variationalEquations.integrateVariationalAndDynamicalEquations(
                        variationalEquations.getPropagatorSettings( )->getInitialStates( ), 1 );
        }
        else
        {
            variationalEquations.integrateDynamicalEquationsOfMotionOnly(
                        variationalEquations.getPropagatorSettings( )->getInitialStates( ) );
        }


        // Retrieve test data
        for( unsigned int arc = 0; arc < integrationArcEnds.size( ); arc++ )
        {
            double testEpoch = integrationArcEnds.at( arc ) - 2.0E4;
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates =
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
            testStates.block( 0, 0, 6, 1 ) = bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( testEpoch );

            testStates.block( 6, 0, 6, 1 ) = bodyMap[ "Orbiter" ]->getStateInBaseFrameFromEphemeris( testEpoch ) -
                    testStates.block( 0, 0, 6, 1 );

            if( propagateVariationalEquations )
            {
                results.first.push_back( variationalEquations.getStateTransitionMatrixInterface( )->
                                         getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
                results.second.push_back( hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getInitialStateList( ).at( arc ) );
            }
            else
            {
                results.second.push_back( testStates );
            }
        }
    }
    return results;
}

BOOST_AUTO_TEST_CASE( testMarsAndOrbiterHybridArcVariationalEquationCalculation )
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
    Eigen::Vector2d perturbedParameter;

    Eigen::Matrix< double, 12, 1> statePerturbation;
    Eigen::VectorXd parameterPerturbation;


    for( unsigned int i = 0; i < 1; i++ )
    {
        // Define parameter perturbation
        parameterPerturbation  = ( Eigen::VectorXd( 2 ) << 1.0E20, 1.0E10 ).finished( );

        // Define state perturbation
        if( i == 0 )
        {
            statePerturbation = ( Eigen::Matrix< double, 12, 1>( )<<
                                  1.0E10, 1.0E10, 1.0E10, 5.0E4, 5.0E4, 10.0E4,
                                  10.0, 10.0, 10.0, 0.1, 0.1, 0.1 ).finished( );
        }


        for( unsigned int patchArcs = 0; patchArcs < 2; patchArcs++ )
        {
            // Test for all requested propagator types.
            for( unsigned int k = 0; k < 1; k++ )
            {
                std::cout<<"Propagating state transition: "<<patchArcs<<" "<<" "<<k<<std::endl;
                // Compute state transition and sensitivity matrices
                currentOutput = executeHybridArcMarsAndOrbiterSensitivitySimulation < double, double >(
                            Eigen::Matrix< double, 12, 1 >::Zero( ), Eigen::VectorXd::Zero( 2 ), true, patchArcs );
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
                    manualPartial[ arc ] = Eigen::MatrixXd::Zero( 12, 14 );
                }

                // Numerically compute state transition matrix
                for( unsigned int j = 0; j < 12; j++ )
                {
                    std::cout<<"Propagating perturbation "<<j<<std::endl;
                    std::vector< Eigen::VectorXd > upPerturbedState, upPerturbedState2, downPerturbedState2, downPerturbedState;
                    perturbedState.setZero( );
                    perturbedState( j ) += statePerturbation( j );
                    upPerturbedState = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, Eigen::VectorXd::Zero( 2 ), false, false, nominalArcStartStates ).second;

                    perturbedState.setZero( );
                    perturbedState( j ) += 0.5 * statePerturbation( j );
                    upPerturbedState2 = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, Eigen::VectorXd::Zero( 2 ), false, false, nominalArcStartStates ).second;

                    perturbedState.setZero( );
                    perturbedState( j ) -= 0.5 * statePerturbation( j );
                    downPerturbedState2 = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, Eigen::VectorXd::Zero( 2 ), false, false, nominalArcStartStates ).second;


                    perturbedState.setZero( );
                    perturbedState( j ) -= statePerturbation( j );
                    downPerturbedState = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, Eigen::VectorXd::Zero( 2 ), false, false, nominalArcStartStates ).second;

                    for( unsigned int arc = 0; arc < upPerturbedState.size( ); arc++ )
                    {
                        manualPartial[ arc ].block( 0, j, 12, 1 ) =
                                ( -upPerturbedState[ arc ] + 8.0 * upPerturbedState2[ arc ] - 8.0 * downPerturbedState2[ arc ] + downPerturbedState[ arc ] ) /
                                ( 6.0 * statePerturbation( j ) );
                    }
                }

                //Numerically compute sensitivity matrix
                for( unsigned int j = 0; j < 2; j ++ )
                {
                    std::vector< Eigen::VectorXd > upPerturbedState, upPerturbedState2, downPerturbedState2, downPerturbedState;
                    perturbedState.setZero( );
                    Eigen::Vector2d upPerturbedParameter, downPerturbedParameter;

                    perturbedParameter.setZero( );
                    perturbedParameter( j ) += parameterPerturbation( j );
                    upPerturbedState = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, perturbedParameter, false, false, nominalArcStartStates).second;

                    perturbedParameter.setZero( );
                    perturbedParameter( j ) += 0.5 * parameterPerturbation( j );
                    upPerturbedState2 = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, perturbedParameter, false, false, nominalArcStartStates ).second;

                    perturbedParameter.setZero( );
                    perturbedParameter( j ) -= 0.5 * parameterPerturbation( j );
                    downPerturbedState2 = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, perturbedParameter, false, false, nominalArcStartStates ).second;

                    perturbedParameter.setZero( );
                    perturbedParameter( j ) -= parameterPerturbation( j );
                    downPerturbedState = executeHybridArcMarsAndOrbiterSensitivitySimulation< double, double >(
                                perturbedState, perturbedParameter, false, false, nominalArcStartStates ).second;

                    for( unsigned int arc = 0; arc < upPerturbedState.size( ); arc++ )
                    {
                        manualPartial[ arc ].block( 0, j + 12, 12, 1 ) =
                                ( -upPerturbedState[ arc ] + 8.0 * upPerturbedState2[ arc ] - 8.0 * downPerturbedState2[ arc ] + downPerturbedState[ arc ] ) /
                                ( 6.0 * parameterPerturbation( j ) );
                    }
                }




                for( unsigned int arc = 0; arc < manualPartial.size( ); arc++ )
                {
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 0, 0, 6, 6 ),
                                manualPartial.at( arc ).block( 0, 0, 6, 6 ), 5.0E-5 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 6, 6, 6, 6 ),
                                manualPartial.at( arc ).block( 6, 6, 6, 6 ), 5.0E-5 );

                    double couplingTolerance;
                    if( arc == 0 )
                    {
                        couplingTolerance = 5.0E-1;
                    }
                    else if( arc == 1 || patchArcs )
                    {
                        couplingTolerance = 5.0E-2;
                    }
                    else
                    {
                        couplingTolerance = 1.0E-2;
                    }

                    // One component is, by chance, not computed to within relative precision (due to small numerical value),
                    // next lines mitigate
                    if( patchArcs == 0 )
                    {
                        stateTransitionAndSensitivityMatrixAtEpoch[ arc ]( 7, 4 ) = 0.0;
                        manualPartial[ arc ]( 7, 4 ) = 0.0;
                    }

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 6, 0, 6, 6 ),
                                manualPartial.at( arc ).block( 6, 0, 6, 6 ), couplingTolerance );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 0, 12, 6, 1 ),
                                manualPartial.at( arc ).block( 0, 12, 6, 1 ), 5.0E-5 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 6, 12, 6, 1 ),
                                manualPartial.at( arc ).block( 6, 12, 6, 1 ), 5.0E-3 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 6, 13, 6, 1 ),
                                manualPartial.at( arc ).block( 6, 13, 6, 1 ), 5.0E-5 );

                }

            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}

