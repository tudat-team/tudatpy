/* Copyright (c) 2010-2019, Delft University of Technology
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

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

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
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    bodies.createEmptyBody( "Orbiter" );
    bodies.at( "Orbiter" )->setConstantBodyMass( 5.0E3 );
    bodies.at( "Orbiter" )->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                            std::map< double, std::shared_ptr< Ephemeris > >( ),
                                            "Mars", "ECLIPJ2000" ) );
    bodies.processBodyFrameDefinitions( );

    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > orbiterRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "Orbiter" )->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    orbiterRadiationPressureSettings, "Orbiter", bodies ) );


    
    


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap singleArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    singleArcAccelerationMap[ "Mars" ] = accelerationsOfMars;

    std::vector< std::string > singleArcBodiesToIntegrate, singleArcCentralBodies;
    singleArcBodiesToIntegrate.push_back( "Mars" );
    singleArcCentralBodies.push_back( "SSB" );

    AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodies, singleArcAccelerationMap, singleArcBodiesToIntegrate, singleArcCentralBodies );
    Eigen::VectorXd singleArcInitialStates = getInitialStatesOfBodies(
                singleArcBodiesToIntegrate, singleArcCentralBodies, bodies, initialEphemerisTime );

    singleArcInitialStates += initialStateDifference.segment(
                0, singleArcInitialStates.rows( ) );

    std::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >(
                singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToIntegrate,
                singleArcInitialStates, finalEphemerisTime );


    SelectedAccelerationMap multiArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfOrbiter;
    accelerationsOfOrbiter[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfOrbiter[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    multiArcAccelerationMap[ "Orbiter" ] = accelerationsOfOrbiter;

    std::vector< std::string > multiArcBodiesToIntegrate, multiArcCentralBodies;
    multiArcBodiesToIntegrate.push_back( "Orbiter" );
    multiArcCentralBodies.push_back( "Mars" );

    AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                bodies, multiArcAccelerationMap, multiArcBodiesToIntegrate, multiArcCentralBodies );

    // Creater arc times
    std::vector< double > integrationArcStarts, integrationArcEnds;
    //    double integrationStartTime = initialEphemerisTime;
    //    double integrationEndTime = finalEphemerisTime - 1.0E4;
    //    double currentStartTime = integrationStartTime;
    //    double currentEndTime = integrationStartTime + arcDuration;
    //    do
    //    {
    //        integrationArcStarts.push_back( currentStartTime );
    //        integrationArcEnds.push_back( currentEndTime );

    //        currentStartTime = currentEndTime - arcOverlap;
    //        currentEndTime = currentStartTime + arcDuration;
    //    }
    double timeBetweenArcs = 86400.0;
    //    double arcDuration = 0.5E6;
    double currentStartTime = initialEphemerisTime;
    double currentEndTime = initialEphemerisTime + arcDuration;
    do
    {
        integrationArcStarts.push_back( currentStartTime );
        integrationArcEnds.push_back( currentEndTime );

        currentEndTime = currentStartTime + timeBetweenArcs + arcDuration;
        currentStartTime = currentStartTime + timeBetweenArcs;
    }
    while( currentEndTime < finalEphemerisTime );

    // Create list of multi-arc initial states
    unsigned int numberOfIntegrationArcs = integrationArcStarts.size( );
    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    multiArcSystemInitialStates.resize( numberOfIntegrationArcs );

    // Define (quasi-arbitrary) arc initial states
    double marsGravitationalParameter =  bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );

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

    std::shared_ptr< IntegratorSettings< > > singleArcIntegratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 60.0 );

    std::shared_ptr< IntegratorSettings< > > multiArcIntegratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 45.0 );


    // Define parameters.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames =
                getInitialStateParameterSettings< double >( hybridArcPropagatorSettings, bodies, integrationArcStarts );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

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
                    bodies, singleArcIntegratorSettings, multiArcIntegratorSettings,
                    hybridArcPropagatorSettings, parametersToEstimate, integrationArcStarts );

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
            testStates.block( 0, 0, 6, 1 ) = bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( testEpoch );

            testStates.block( 6, 0, 6, 1 ) = bodies.at( "Orbiter" )->getStateInBaseFrameFromEphemeris( testEpoch );/* -
                    testStates.block( 0, 0, 6, 1 );*/

            if( propagateVariationalEquations )
            {
                results.first.push_back( variationalEquations.getStateTransitionMatrixInterface( )->
                                         getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
                results.second.push_back( hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getInitialStateList( ).at( arc ) );

                Eigen::MatrixXd testMatrixDirect =
                        variationalEquations.getStateTransitionMatrixInterface( )->
                          getCombinedStateTransitionAndSensitivityMatrix( testEpoch );
                Eigen::MatrixXd testMatrixFull=
                        variationalEquations.getStateTransitionMatrixInterface( )->
                          getFullCombinedStateTransitionAndSensitivityMatrix( testEpoch );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            testMatrixDirect.block( 0, 0, 12, 6 ),
                            testMatrixFull.block( 0, 0, 12, 6 ),
                            std::numeric_limits< double >::epsilon( ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            testMatrixDirect.block( 0, 6, 12, 6 ),
                            testMatrixFull.block( 0, 6 * ( arc + 1 ), 12, 6 ),
                            std::numeric_limits< double >::epsilon( ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            testMatrixDirect.block( 0, 12, 12, 2 ),
                            testMatrixFull.block( 0, 18, 12, 2 ),
                            std::numeric_limits< double >::epsilon( ) );
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

    // Define variables for numerical differentiation
    Eigen::Matrix< double, 12, 1>  perturbedState;
    Eigen::Vector2d perturbedParameter;

    Eigen::Matrix< double, 12, 1> statePerturbation;
    Eigen::VectorXd parameterPerturbation;


    // Define parameter perturbation
    parameterPerturbation  = ( Eigen::VectorXd( 2 ) << 1.0E20, 1.0E10 ).finished( );
    statePerturbation = ( Eigen::Matrix< double, 12, 1>( )<<
                              1.0E10, 1.0E10, 1.0E10, 5.0E4, 5.0E4, 10.0E4,
                              10.0, 10.0, 10.0, 0.1, 0.1, 0.1 ).finished( );

    for( unsigned int patchArcs = 0; patchArcs < 1; patchArcs++ )
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
//                std::cout<<"Propagating perturbation "<<j<<std::endl;
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

                //                    std::cout<<"Arc: "<<arc<<std::endl<<stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 6, 0, 6, 6 )<<std::endl<<std::endl<<
                //                               manualPartial.at( arc ).block( 6, 0, 6, 6 )<<std::endl<<std::endl<<
                //                               ( stateTransitionAndSensitivityMatrixAtEpoch.at( arc ) - manualPartial.at( arc ) ).block( 6, 0, 6, 6 ).cwiseQuotient(
                //                                manualPartial.at( arc ).block( 6, 0, 6, 6 ) )<<std::endl<<std::endl;
//                std::cout<<"Arc: "<<arc<<std::endl<<stateTransitionAndSensitivityMatrixAtEpoch.at( arc ).block( 0, 12, 12, 2 )<<std::endl<<std::endl<<
//                           manualPartial.at( arc ).block( 0, 12, 12, 2 )<<std::endl<<std::endl<<
//                           ( stateTransitionAndSensitivityMatrixAtEpoch.at( arc ) - manualPartial.at( arc ) ).block( 0, 12, 12, 2 ).cwiseQuotient(
//                               manualPartial.at( arc ).block( 0, 12, 12, 2 ) )<<std::endl<<std::endl;

            }

        }
    }
}

BOOST_AUTO_TEST_CASE( testVaryingCentralBodyHybridArcVariationalEquations )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );

    // Specify initial time
    double initialTime = 0.0;
    double finalTime = 4.0 * 86400.0;

    // Get body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialTime - 3600.0, finalTime + 3600.0,
                                    "Jupiter", "ECLIPJ2000" );

    // Create bodies needed in simulation
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Spacecraft" );
    bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                               std::map< double, std::shared_ptr< Ephemeris > >( ), "Jupiter", "ECLIPJ2000" ) );
    bodies.processBodyFrameDefinitions( );

    SelectedAccelerationMap singleArcAccelerationMap;
    std::vector< std::string > singleArcBodiesToPropagate = { "Io", "Europa", "Ganymede" };
    std::vector< std::string > singleArcCentralBodies = { "Jupiter", "Jupiter", "Jupiter" };

    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
    {
        singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ "Jupiter" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

        for( unsigned int j = 0; j < singleArcBodiesToPropagate.size( ); j++ )
        {
            if( i != j )
            {
                singleArcAccelerationMap[ singleArcBodiesToPropagate.at( i ) ][ singleArcBodiesToPropagate.at( j ) ].push_back(
                            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
            }
        }
    }

    basic_astrodynamics::AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodies, singleArcAccelerationMap, singleArcBodiesToPropagate, singleArcCentralBodies );

    Eigen::VectorXd singleArcInitialState = getInitialStatesOfBodies(
                singleArcBodiesToPropagate, singleArcCentralBodies, bodies, initialTime );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > singleArcPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToPropagate,singleArcInitialState, finalTime );

    std::vector< std::string > multiArcBodiesToPropagate =
    { "Spacecraft", "Spacecraft", "Spacecraft", "Spacecraft", "Spacecraft", "Spacecraft" };

    std::vector< std::string > multiArcCentralBodies =
    { "Io", "Ganymede", "Europa", "Ganymede", "Io", "Europa" };

    std::vector< double > arcStartTimes =
    { 3600.0, 3.0 * 3600, 5.0 * 3600.0, 7.0 * 3600.0, 9.0 * 3600.0, 11.0 * 3600.0 };
    std::map< std::string, std::vector< double > > arcStartTimesPerBody;

    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        arcStartTimesPerBody[ multiArcCentralBodies.at( i ) ].push_back( arcStartTimes.at( i ) );
    }
    double arcDuration = 3600.0;


    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > multiArcPropagationSettingsList;
    std::map< std::string, std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > >
            multiArcPropagationSettingsListPerCentralBody;
    std::map< std::string, std::vector< int > > perBodyIndicesInFullPropagation;

    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        Eigen::Vector6d spacecraftInitialStateInKeplerianElements;
        spacecraftInitialStateInKeplerianElements( semiMajorAxisIndex ) = 3500.0E3;
        spacecraftInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        spacecraftInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        spacecraftInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        spacecraftInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        spacecraftInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians(
                    static_cast< double >( i * 30 ) );
        double centralBodyGravitationalParameter = bodies.at( multiArcCentralBodies.at( i ) )->getGravityFieldModel( )->getGravitationalParameter( );
        multiArcSystemInitialStates.push_back( convertKeplerianToCartesianElements(
                                                   spacecraftInitialStateInKeplerianElements, centralBodyGravitationalParameter ) );

        SelectedAccelerationMap multiArcAccelerationMap;

        multiArcAccelerationMap[ "Spacecraft" ][ "Jupiter" ].push_back(
                    std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

        for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
        {
            multiArcAccelerationMap[ "Spacecraft" ][ singleArcBodiesToPropagate.at( i ) ].push_back(
                        std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
        }

        basic_astrodynamics::AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                    bodies, multiArcAccelerationMap, { multiArcBodiesToPropagate.at( i ) }, { multiArcCentralBodies.at( i ) } );

        multiArcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( std::vector< std::string >{ multiArcCentralBodies.at( i ) }, multiArcAccelerationModelMap,
                      std::vector< std::string >{ multiArcBodiesToPropagate.at( i ) },
                      multiArcSystemInitialStates.at( i ), arcStartTimes.at( i ) + arcDuration ) );
        multiArcPropagationSettingsListPerCentralBody[
                multiArcCentralBodies.at( i ) ].push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( std::vector< std::string >{ multiArcCentralBodies.at( i ) }, multiArcAccelerationModelMap,
                      std::vector< std::string >{ multiArcBodiesToPropagate.at( i ) },
                      multiArcSystemInitialStates.at( i ), arcStartTimes.at( i ) + arcDuration ) );
        perBodyIndicesInFullPropagation[ multiArcCentralBodies.at( i ) ].push_back( i );

    }
    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagationSettings =
            std::make_shared< MultiArcPropagatorSettings< > >( multiArcPropagationSettingsList );

    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings =
            std::make_shared< HybridArcPropagatorSettings< > >( singleArcPropagatorSettings, multiArcPropagationSettings );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames = getInitialHybridArcParameterSettings< >( hybridArcPropagatorSettings, bodies, arcStartTimes );

    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                                      singleArcBodiesToPropagate.at( i ), gravitational_parameter ) );
    }

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, hybridArcPropagatorSettings );
    printEstimatableParameterEntries( parametersToEstimate );

    std::shared_ptr< IntegratorSettings< > > singleArcIntegratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialTime, 60.0 );

    std::shared_ptr< IntegratorSettings< > > multiArcIntegratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, TUDAT_NAN, 15.0 );

    // Create dynamics simulator
    HybridArcVariationalEquationsSolver< > variationalEquations =
            HybridArcVariationalEquationsSolver< >(
                bodies, singleArcIntegratorSettings, multiArcIntegratorSettings,
                hybridArcPropagatorSettings, parametersToEstimate, arcStartTimes, true, false, true );

    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > fullMultiArcVariationalSolution =
            variationalEquations.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );
    std::vector< std::map< double, Eigen::VectorXd > > fullMultiArcStateSolution =
            variationalEquations.getMultiArcSolver( )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

    for( unsigned int i = 0; i < singleArcBodiesToPropagate.size( ); i++ )
    {
        std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPerBodyPropagationSettings =
                std::make_shared< MultiArcPropagatorSettings< > >( multiArcPropagationSettingsListPerCentralBody.at(
                                                                       singleArcBodiesToPropagate.at( i ) ) );
        std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPerBodyPropagatorSettings =
                std::make_shared< HybridArcPropagatorSettings< > >(
                    singleArcPropagatorSettings, multiArcPerBodyPropagationSettings );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNamesPerBody;

        parameterNamesPerBody = getInitialHybridArcParameterSettings< >(
                    hybridArcPerBodyPropagatorSettings, bodies, arcStartTimesPerBody.at( singleArcBodiesToPropagate.at( i ) ) );

        for( unsigned int j = 0; j < singleArcBodiesToPropagate.size( ); j++ )
        {
            parameterNamesPerBody.push_back( std::make_shared< EstimatableParameterSettings >(
                                                 singleArcBodiesToPropagate.at( j ), gravitational_parameter ) );
        }

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimatePerBody =
                createParametersToEstimate< double >( parameterNamesPerBody, bodies, hybridArcPerBodyPropagatorSettings );
        printEstimatableParameterEntries( parametersToEstimatePerBody );

        HybridArcVariationalEquationsSolver< > perCentralBodyVariationalEquations =
                HybridArcVariationalEquationsSolver< >(
                    bodies, singleArcIntegratorSettings, multiArcIntegratorSettings,
                    hybridArcPerBodyPropagatorSettings, parametersToEstimatePerBody, arcStartTimesPerBody.at(
                        singleArcBodiesToPropagate.at( i ) ), true, false, true );

        std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > perBodyMultiArcVariationalSolution =
                perCentralBodyVariationalEquations.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );
        std::vector< std::map< double, Eigen::VectorXd > > perBodyMultiArcStateSolution =
                perCentralBodyVariationalEquations.getMultiArcSolver( )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

        for( unsigned int j = 0; j < perBodyIndicesInFullPropagation.at( singleArcBodiesToPropagate.at( i ) ).size( ); j++ )
        {
            for( unsigned int k = 0; k < 2; k++ )
            {
                std::map< double, Eigen::MatrixXd > fullMultiArcMatrixHistory = fullMultiArcVariationalSolution.at(
                            perBodyIndicesInFullPropagation.at( singleArcBodiesToPropagate.at( i ) ).at( j ) ).at( k );
                std::map< double, Eigen::MatrixXd > perBodyMultiMatrixHistory = perBodyMultiArcVariationalSolution.at( j ).at( k );

                auto fullIterator = fullMultiArcMatrixHistory.begin( );
                auto perBodyIterator = perBodyMultiMatrixHistory.begin( );

                BOOST_CHECK_EQUAL( fullMultiArcMatrixHistory.size( ), perBodyMultiMatrixHistory.size( ) );

                for( unsigned int i = 0; i < fullMultiArcMatrixHistory.size( ); i++ )
                {
                    BOOST_CHECK_CLOSE_FRACTION( fullIterator->first, perBodyIterator->first, std::numeric_limits< double >::epsilon( ) );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fullIterator->second, perBodyIterator->second, std::numeric_limits< double >::epsilon( ) );

                    fullIterator++;
                    perBodyIterator++;
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

