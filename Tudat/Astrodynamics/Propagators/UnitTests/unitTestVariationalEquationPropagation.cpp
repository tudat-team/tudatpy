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
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"

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

BOOST_AUTO_TEST_SUITE( test_variational_equation_calculation )

template< typename TimeType = double , typename StateScalarType  = double >
        std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
executeEarthMoonSimulation(
        const std::vector< std::string > centralBodies,
        const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference =
        Eigen::Matrix< StateScalarType, 12, 1 >::Zero( ),
        const int propagationType = 0,
        const Eigen::Vector3d parameterPerturbation = Eigen::Vector3d::Zero( ),
        const bool propagateVariationalEquations = 1 )
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
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Earth" ] = accelerationsOfEarth;

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    accelerationsOfMoon[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
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
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime ), 1800.0 );


    // Set initial states of bodies to integrate.
    TimeType initialIntegrationTime = initialEphemerisTime;

    // Set (perturbed) initial state.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialTranslationalState;
    initialTranslationalState = getInitialStatesOfBodies< TimeType, StateScalarType >(
                bodiesToIntegrate, centralBodies, bodyMap, initialIntegrationTime );
    initialTranslationalState += initialStateDifference;

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings;
    TranslationalPropagatorType propagatorType;
    if( propagationType == 0 )
    {
        propagatorType = cowell;
    }
    else if( propagationType == 1 )
    {
        propagatorType = encke;
    }
    propagatorSettings =  std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialTranslationalState,
              TimeType( finalEphemerisTime ), propagatorType );

    // Define parameters.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames.push_back(
                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        "Moon", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                            "Moon", centralBodies[ 0 ], bodyMap, TimeType( initialEphemerisTime ) ) +
                    initialStateDifference.segment( 0, 6 ),
                    centralBodies[ 0 ] ) );
        parameterNames.push_back(
                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        "Earth", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                            "Earth", centralBodies[ 1 ], bodyMap, TimeType( initialEphemerisTime ) ) +
                    initialStateDifference.segment( 6, 6 ),
                    centralBodies[ 1 ] ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );

    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Perturb parameters.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    parameterVector.block( 12, 0, 3, 1 ) += parameterPerturbation;
    parametersToEstimate->resetParameterValues( parameterVector );

    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;

    {
        // Create dynamics simulator
        SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
                SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
                    1, std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), 1, 0 );

        // Propagate requested equations.
        if( propagateVariationalEquations )
        {
            dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
        }
        else
        {
            dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
        }

        // Retrieve test data
        double testEpoch = 1.4E7;
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
            results.first.push_back( dynamicsSimulator.getStateTransitionMatrixInterface( )->
                                     getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
        }
        results.second.push_back( testStates );
    }
    return results;
}

//! Test the state transition and sensitivity matrix computation against their numerical propagation.
/*!
 *  Test the state transition and sensitivity matrix computation against their numerical propagation. This unit test
 *  propagates the variational equations for the Earth and Moon, using both a barycentric origin and hierarchical origin.
 *  For the hierarchical origin, both an Encke and Cowell propagator are used. The results are compared against
 *  results obtained from numerical differentiation (first-order central difference). *
 */
BOOST_AUTO_TEST_CASE( testEarthMoonVariationalEquationCalculation )
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


    for( unsigned int i = 0; i < centralBodiesSet.size( ); i++ )
    {
        // Define parameter perturbation
        parameterPerturbation  = ( Eigen::Vector3d( ) << 1.0E10, 1.0E10, 1.0E14 ).finished( );

        // Define state perturbation
        if( i == 0 )
        {
            statePerturbation = ( Eigen::Matrix< double, 12, 1>( ) <<
                                  100000.0, 100000.0, 100000.0, 0.1, 0.1, 0.1,
                                  100000.0, 100000.0, 100000.0, 0.1, 0.1, 0.1 ).finished( );
        }
        else if( i == 1 )
        {
            statePerturbation = ( Eigen::Matrix< double, 12, 1>( ) <<
                                  100000.0, 100000.0, 100000.0, 0.1, 0.1, 0.1,
                                  100000.0, 100000.0, 10000000.0, 0.1, 0.1, 10.0 ).finished( );
        }

        // Only perform Encke propagation is origin is not SSB
        unsigned int maximumPropagatorType = 1;
        if( i == 1 )
        {
            maximumPropagatorType = 2;
        }

        // Test for all requested propagator types.
        for( unsigned int k = 0; k < maximumPropagatorType; k++ )
        {
            // Compute state transition and sensitivity matrices
            currentOutput = executeEarthMoonSimulation< double, double >(
                        centralBodiesSet[ i ], Eigen::Matrix< double, 12, 1 >::Zero( ), k );
            Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.at( 0 );
            Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 12, 15 );

            // Numerically compute state transition matrix
            for( unsigned int j = 0; j < 12; j++ )
            {
                Eigen::VectorXd upPerturbedState, downPerturbedState;
                perturbedState.setZero( );
                perturbedState( j ) += statePerturbation( j );
                upPerturbedState = executeEarthMoonSimulation< double, double >(
                            centralBodiesSet[ i ], perturbedState, k, Eigen::Vector3d::Zero( ), 0 ).second.at( 0 );

                perturbedState.setZero( );
                perturbedState( j ) -= statePerturbation( j );
                downPerturbedState = executeEarthMoonSimulation< double, double >(
                            centralBodiesSet[ i ], perturbedState, k, Eigen::Vector3d::Zero( ), 0 ).second.at( 0 );
                manualPartial.block( 0, j, 12, 1 ) =
                        ( upPerturbedState - downPerturbedState ) / ( 2.0 * statePerturbation( j ) );
            }

            // Numerically compute sensitivity matrix
            for( unsigned int j = 0; j < 3; j ++ )
            {
                Eigen::VectorXd upPerturbedState, downPerturbedState;
                perturbedState.setZero( );
                Eigen::Vector3d upPerturbedParameter, downPerturbedParameter;
                perturbedParameter.setZero( );
                perturbedParameter( j ) += parameterPerturbation( j );
                upPerturbedState = executeEarthMoonSimulation< double, double >(
                            centralBodiesSet[ i ], perturbedState, k, perturbedParameter, 0 ).second.at( 0 );

                perturbedParameter.setZero( );
                perturbedParameter( j ) -= parameterPerturbation( j );
                downPerturbedState = executeEarthMoonSimulation< double, double >(
                            centralBodiesSet[ i ], perturbedState, k, perturbedParameter, 0 ).second.at( 0 );
                manualPartial.block( 0, j + 12, 12, 1 ) =
                        ( upPerturbedState - downPerturbedState ) / ( 2.0 * parameterPerturbation( j ) );
            }

            // Check results
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 12, 15 ), manualPartial, 2.0E-4 );

        }

    }
}

template< typename TimeType = double , typename StateScalarType  = double >
std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
executeOrbiterSimulation(
        const Eigen::Matrix< StateScalarType, 6, 1 > initialStateDifference =
        Eigen::Matrix< StateScalarType, 6, 1 >::Zero( ),
        const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 10 ),
        const bool propagateVariationalEquations = 1 )
{
    int numberOfParametersToEstimate = 10;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = initialEphemerisTime + 4.0 * 3600.0;

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames );
    NamedBodyMap bodyMap = createBodies( bodySettings );
    bodyMap[ "Vehicle" ] = std::make_shared< Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.1, -0.4 ).finished( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Vehicle" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Vehicle", bodyMap ) );

    bodyMap[ "Vehicle" ]->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                            < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );

    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::cannon_ball_radiation_pressure ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;


    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToIntegrate, centralBodies );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime ), 5.0 );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< StateScalarType, 6, 1 > initialTranslationalState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );
    initialTranslationalState += initialStateDifference;

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialTranslationalState,
              TimeType( finalEphemerisTime ), cowell );

    // Define parameters.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames.push_back(
                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        "Vehicle", initialTranslationalState, "Earth" ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      3, 0, 3, 3, "Earth", spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      3, 1, 3, 3, "Earth", spherical_harmonics_sine_coefficient_block ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Perturb parameters.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    parameterVector.block( 6, 0, numberOfParametersToEstimate, 1 ) += parameterPerturbation;
    parametersToEstimate->resetParameterValues( parameterVector );


    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;

    {
        // Create dynamics simulator
        SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
                SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
                    1, std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), 0, 0 );

        // Propagate requested equations.
        if( propagateVariationalEquations )
        {
            dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
        }
        else
        {
            dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
        }

        // Retrieve test data
        double testEpoch = initialEphemerisTime + 3.5 * 3600.0;
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 12 );
        testStates.block( 0, 0, 6, 1 ) = bodyMap[ "Vehicle" ]->getEphemeris( )->getCartesianState( testEpoch );

        if( propagateVariationalEquations )
        {
            results.first.push_back( dynamicsSimulator.getStateTransitionMatrixInterface( )->
                                     getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );
        }
        results.second.push_back( testStates );
    }
    return results;
}

//! Test the state transition and sensitivity matrix computation against their numerical propagation for Earth orbiter.
/*!
 *  Test the state transition and sensitivity matrix computation against their numerical propagation for Earth orbiter.
 *  This unit test propagates the variational equations for a low Earth orbiter, using the Enckle propagator.
 *  The results are compared against results obtained from numerical differentiation (first-order central difference).
 *  The estimated parameters consist of a mass of a third-body, radiation pressure coefficient and drag coefficient of the
 *  vehicle and spherical harmonic coefficients of the Earth
 */
BOOST_AUTO_TEST_CASE( testEarthOrbiterVariationalEquationCalculation )
{
    std::pair< std::vector< Eigen::MatrixXd >, std::vector< Eigen::VectorXd > > currentOutput;

    // Define variables for numerical differentiation
    Eigen::Matrix< double, 6, 1 > perturbedState;
    Eigen::Matrix< double, 6, 1 > statePerturbation;

    // Define parameter perturbation
    int numberOfParametersToEstimate = 10;
    double sphericalHarmonicsPerturbation = 1.0E-6;
    Eigen::Matrix< double, 10, 1 > perturbedParameter;
    Eigen::Matrix< double, 10, 1 > parameterPerturbation  =
            ( Eigen::Matrix< double, 10, 1 >( ) << 10.0, 1.0, 1.0E12,
              sphericalHarmonicsPerturbation, sphericalHarmonicsPerturbation,
              sphericalHarmonicsPerturbation, sphericalHarmonicsPerturbation,
              sphericalHarmonicsPerturbation, sphericalHarmonicsPerturbation, sphericalHarmonicsPerturbation ).finished( );

    // Define state perturbation
    statePerturbation = ( Eigen::Matrix< double, 6, 1>( ) <<
                          10.0, 10.0, 10.0, 0.01, 0.01, 0.01 ).finished( );


    // Compute state transition and sensitivity matrices
    currentOutput = executeOrbiterSimulation< double, double >(
                Eigen::Matrix< double, 6, 1 >::Zero( ) );
    Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.at( 0 );
    Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 6, 6 + numberOfParametersToEstimate );

    // Numerically compute state transition matrix
    for( unsigned int j = 0; j < 6; j++ )
    {
        Eigen::VectorXd upPerturbedState, downPerturbedState;
        perturbedState.setZero( );
        perturbedState( j ) += statePerturbation( j );
        upPerturbedState = executeOrbiterSimulation< double, double >(
                    perturbedState, Eigen::Matrix< double, 10, 1 >::Zero( ), 0 ).second.at( 0 );

        perturbedState.setZero( );
        perturbedState( j ) -= statePerturbation( j );
        downPerturbedState = executeOrbiterSimulation< double, double >(
                    perturbedState, Eigen::Matrix< double, 10, 1 >::Zero( ), 0 ).second.at( 0 );

        manualPartial.block( 0, j, 6, 1 ) =
                ( upPerturbedState.segment( 0, 6 ) - downPerturbedState.segment( 0, 6 ) ) / ( 2.0 * statePerturbation( j ) );
    }

    // Numerically compute sensitivity matrix
    for( int j = 0; j < numberOfParametersToEstimate; j ++ )
    {
        Eigen::VectorXd upPerturbedState, downPerturbedState;
        perturbedState.setZero( );
        Eigen::Matrix< double, 10, 1 > upPerturbedParameter, downPerturbedParameter;
        perturbedParameter.setZero( );
        perturbedParameter( j ) += parameterPerturbation( j );

        upPerturbedState = executeOrbiterSimulation< double, double >(
                    perturbedState, perturbedParameter ).second.at( 0 );

        perturbedParameter.setZero( );
        perturbedParameter( j ) -= parameterPerturbation( j );
        downPerturbedState = executeOrbiterSimulation< double, double >(
                    perturbedState, perturbedParameter ).second.at( 0 );

        manualPartial.block( 0, j + 6, 6, 1 ) =
                ( upPerturbedState.segment( 0, 6 ) - downPerturbedState.segment( 0, 6 ) ) / ( 2.0 * parameterPerturbation( j ) );
    }

    // Check results
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                stateTransitionAndSensitivityMatrixAtEpoch, manualPartial, 5.0E-5 );

}

template< typename TimeType = double , typename StateScalarType  = double >
        std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
executePhobosRotationSimulation(
        const Eigen::Matrix< StateScalarType, 13, 1 > initialStateDifference,
        Eigen::Matrix< StateScalarType, 13, 1 >& appliedStateDifference,
        const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 8 ),
        const bool propagateVariationalEquations = 1 )
{
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 86400.0;
    int numberOfParametersToEstimate = 8;

    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = std::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >(
                                         [ = ]( ){ return Eigen::Vector6d::Zero( ); } ) );
    bodyMap[ "Mars" ]->setRotationalEphemeris(
                simulation_setup::createRotationModel( simulation_setup::getDefaultRotationModelSettings(
                                                           "Mars", initialEphemerisTime, finalEphemerisTime ), "Mars" ) );
    bodyMap[ "Mars" ]->setGravityFieldModel(
                simulation_setup::createGravityFieldModel(
                    simulation_setup::getDefaultGravityFieldSettings(
                        "Mars", initialEphemerisTime, finalEphemerisTime ), "Mars", bodyMap ) );
    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );

    bodyMap[ "Phobos" ] = std::make_shared< Body >( );

    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;


    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );
    bodyMap[ "Phobos" ]->setBodyInertiaTensor(
                phobosInertiaTensor, ( 0.3615 + 0.4265 + 0.5024 ) / 3.0 );

    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
            phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed",
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true ) ) );

    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;
    phobosKeplerElements( 2 ) = 0.1;

    bodyMap[ "Phobos" ]->setEphemeris(
                tudat::ephemerides::getTabulatedEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                                                               phobosKeplerElements, 0.0, marsGravitationalParameter ),
                                                           initialEphemerisTime - 3600.0,
                                                           finalEphemerisTime + 3600, 60.0 ) );

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::AngleAxisd( 1.0E-0, Eigen::Vector3d::UnitZ( ) ) *
                                                                  Eigen::AngleAxisd( 2.0E-0, Eigen::Vector3d::UnitX( ) ) *
                                                                  Eigen::AngleAxisd( -0.5E-0, Eigen::Vector3d::UnitZ( ) ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );
    unitRotationState( 4 ) = 1.0E-4;
    unitRotationState( 5 ) = -1.0E-4;
    unitRotationState( 6 ) = std::sqrt( marsGravitationalParameter /
                                        std::pow( phobosSemiMajorAxis, 3.0 ) ) + 1.0E-4;


    Eigen::Matrix< double, 7, 1 > originalRotationState = unitRotationState;
    Eigen::Matrix< double, 7, 1 > stateDifferenceToAdd = initialStateDifference.segment( 6, 7 );

    unitRotationState += stateDifferenceToAdd;
    unitRotationState( 0 ) = originalRotationState( 0 ) / std::fabs( originalRotationState( 0 ) ) *
            std::sqrt( 1.0 - std::pow( unitRotationState.segment( 1, 3 ).norm( ), 2.0 ) );

    appliedStateDifference.segment( 0, 6 ) = initialStateDifference.segment( 0, 6 );
    appliedStateDifference.segment( 6, 7 ) = unitRotationState - originalRotationState ;
    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodyMap[ "Phobos" ]->setRotationalEphemeris( std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

    // Create empty bodies, phobos and mars.
    std::shared_ptr< Body > phobos = bodyMap.at( "Phobos" );
    std::shared_ptr< Body > mars = bodyMap.at( "Mars" );
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
//    accelerationMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
   accelerationMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );

    std::vector< std::string > translationalBodiesToIntegrate;
    std::vector< std::string > translationalCentralBodies;

    translationalBodiesToIntegrate.push_back( "Phobos" );
    translationalCentralBodies.push_back( "Mars" );


    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, translationalBodiesToIntegrate, translationalCentralBodies );


    SelectedTorqueMap torqueMap;
    torqueMap[ "Phobos" ][ "Mars" ].push_back(
                std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );

    // Define propagator settings.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Phobos" );

    // Create torque models
    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodyMap, torqueMap, bodiesToIntegrate );


    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToIntegrate, unitRotationState, std::make_shared< PropagationTimeTerminationSettings >(
                  finalEphemerisTime ) );

    Eigen::VectorXd initialTranslationalState;
    initialTranslationalState = getInitialStatesOfBodies(
                translationalBodiesToIntegrate, translationalCentralBodies, bodyMap, initialEphemerisTime );

    initialTranslationalState += initialStateDifference.segment( 0, 6 );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >
            ( translationalCentralBodies, accelerationModelMap, translationalBodiesToIntegrate, initialTranslationalState,
              finalEphemerisTime, cowell );


    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
    propagatorSettingsList.push_back( rotationalPropagatorSettings );
    propagatorSettingsList.push_back( translationalPropagatorSettings );

    std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsList,
                std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ) );



    // Create integrator settings
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime ), 15.0 );

    // Define parameters.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    {
        parameterNames.push_back( std::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >(
                                      "Phobos", unitRotationState, "ECLIPJ2000" ) );
        parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                      "Phobos", initialTranslationalState, "Mars" ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      1, 0, 2, 2, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      2, 1, 2, 2, "Phobos", spherical_harmonics_sine_coefficient_block ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    Eigen::MatrixXd constraintStateMultiplier;
    Eigen::VectorXd constraintRightHandSide;
    parametersToEstimate->getConstraints( constraintStateMultiplier, constraintRightHandSide );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( constraintStateMultiplier.block( 0, 0, 1, 4 ) ), ( unitRotationState.segment( 0, 4 ) ).transpose( ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( constraintStateMultiplier.block( 0, 4, 1, 9 ) ), ( Eigen::MatrixXd::Zero( 1, 9 ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( constraintRightHandSide.block( 0, 0, 1, 1 ) ), ( Eigen::MatrixXd::Zero( 1, 1 ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Perturb parameters.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterVector =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    parameterVector.block( 13, 0, numberOfParametersToEstimate, 1 ) += parameterPerturbation;
//    std::cout<<"Parameter perturbation "<<
//               ( parametersToEstimate->template getFullParameterValues< StateScalarType >( ) -
//                 parameterVector ).transpose( )<<std::endl<<
//               parameterPerturbation.transpose( )<<std::endl;;
    parametersToEstimate->resetParameterValues( parameterVector );


    std::pair< std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > results;

    {
        // Create dynamics simulator
        propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > dynamicsSimulator =
                propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType >(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
                    1, std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), 0, 0 );

        // Propagate requested equations.
        if( propagateVariationalEquations )
        {
            dynamicsSimulator.integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
        }
        else
        {
            dynamicsSimulator.integrateDynamicalEquationsOfMotionOnly( propagatorSettings->getInitialStates( ) );
        }


        //        tudat::input_output::writeDataMapToTextFile(
        //                    dynamicsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ),
        //                    "rotPropTest.dat" );



        // Retrieve test data
        double testEpoch = initialEphemerisTime + 15.0 * 3600.0;
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > testStates =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 13 );
        testStates.segment( 6, 7 ) = bodyMap[ "Phobos" ]->getRotationalEphemeris( )->getRotationStateVector( testEpoch );
        testStates.segment( 0, 6 ) = bodyMap[ "Phobos" ]->getEphemeris( )->getCartesianState( testEpoch );

        if( propagateVariationalEquations )
        {
            results.first.push_back( dynamicsSimulator.getStateTransitionMatrixInterface( )->
                                     getCombinedStateTransitionAndSensitivityMatrix( testEpoch ) );

        }
        results.second.push_back( testStates );


    }
    return results;
}

BOOST_AUTO_TEST_CASE( testPhobosRotationVariationalEquationCalculation )
{
    spice_interface::loadStandardSpiceKernels( );

    std::pair< std::vector< Eigen::MatrixXd >, std::vector< Eigen::VectorXd > > currentOutput;

    // Define variables for numerical differentiation
    Eigen::Matrix< double, 13, 1 > perturbedState;
    Eigen::Matrix< double, 13, 1 > statePerturbation;

    // Define parameter perturbation
    int numberOfParametersToEstimate = 8;
    double sphericalHarmonicsPerturbation = 1.0E-4;
    Eigen::Matrix< double, 8, 1 > perturbedParameter;
    Eigen::Matrix< double, 8, 1 > parameterPerturbation;
    parameterPerturbation = Eigen::Matrix< double, 8, 1 >::Constant( sphericalHarmonicsPerturbation );
    parameterPerturbation( 0 ) = 1.0E-4;

    // Compute state transition and sensitivity matrices
    Eigen::Matrix< double, 13, 1 > appliedStateDifference;
    currentOutput = executePhobosRotationSimulation< double, double >(
                Eigen::Matrix< double, 13, 1 >::Zero( ), appliedStateDifference );
    Eigen::MatrixXd stateTransitionAndSensitivityMatrixAtEpoch = currentOutput.first.at( 0 );
    Eigen::VectorXd nominalState = currentOutput.second.at( 0 );

    // Define state perturbation
    statePerturbation = ( Eigen::Matrix< double, 13, 1>( ) <<
                          10.0, 10.0, 10.0, 0.1, 0.01, 0.01,
                          1.0E-5, 1.0E-5, 1.0E-5, 1.0E-5, 1.0E-7, 1.0E-7, 1.0E-7 ).finished( );

    Eigen::MatrixXd manualPartial = Eigen::MatrixXd::Zero( 13, 13 + numberOfParametersToEstimate );

    // Numerically compute state transition matrix
    for( unsigned int test = 0; test < 2; test++ )
    {
        double perturbationMultiplier = ( test == 0 ? 1.0 : 1.0E-3 );
        for( unsigned int j = 7; j < 10; j++ )
        {
            Eigen::Matrix< double, 13, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

            Eigen::VectorXd upPerturbedState, downPerturbedState;
            perturbedState.setZero( );
            perturbedState( j ) += perturbationMultiplier * statePerturbation( j );
            upPerturbedState = executePhobosRotationSimulation< double, double >(
                        perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.at( 0 );

            Eigen::VectorXd stateDifferenceUp = upPerturbedState - nominalState;

            std::cout<<"Test output "<<test<<" "<<j<<"stateDifferenceUp"<<std::endl<<
                       stateDifferenceUp<<std::endl<<std::endl<<
                       "stateTransitionAndSensitivityMatrixAtEpoch"<<std::endl<<
                                              stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
                       "appliedStateDifferenceUp"<<std::endl<<
                                              appliedStateDifferenceUp<<std::endl<<std::endl<<
                       "( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp )"<<std::endl<<
                                              ( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp )<<std::endl<<std::endl;
            if( test == 0 )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( ( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp ).segment( 0, 6 ) ),
                            ( stateDifferenceUp.segment( 0, 6 ) ), 1.0E-3 );
            }
            else
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( ( stateTransitionAndSensitivityMatrixAtEpoch * appliedStateDifferenceUp ).segment( 6, 7 ) ),
                            ( stateDifferenceUp.segment( 6, 7 ) ), 1.0E-5 );
            }
        }
    }

    for( unsigned int j = 0; j < 6; j++ )
    {
        Eigen::Matrix< double, 13, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

        Eigen::VectorXd upPerturbedState, downPerturbedState;
        perturbedState.setZero( );
        perturbedState( j ) += statePerturbation( j );
        upPerturbedState = executePhobosRotationSimulation< double, double >(
                    perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.at( 0 );

        perturbedState.setZero( );
        perturbedState( j ) -= statePerturbation( j );
        downPerturbedState = executePhobosRotationSimulation< double, double >(
                    perturbedState, appliedStateDifferenceDown, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.at( 0 );

        manualPartial.block( 0, j, 13, 1 ) =
                ( upPerturbedState.segment( 0, 13 ) - downPerturbedState.segment( 0, 13 ) ) / ( 2.0 * statePerturbation( j ) );

    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualPartial.block( 0, 0, 13, 6 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 0, 13, 6 ) ), 1.0E-4 );

    for( unsigned int j = 10; j < 13; j++ )
    {
        Eigen::Matrix< double, 13, 1 > appliedStateDifferenceUp, appliedStateDifferenceDown;

        Eigen::VectorXd upPerturbedState, downPerturbedState;
        perturbedState.setZero( );
        perturbedState( j ) += statePerturbation( j );
        upPerturbedState = executePhobosRotationSimulation< double, double >(
                    perturbedState, appliedStateDifferenceUp, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.at( 0 );

        perturbedState.setZero( );
        perturbedState( j ) -= statePerturbation( j );
        downPerturbedState = executePhobosRotationSimulation< double, double >(
                    perturbedState, appliedStateDifferenceDown, Eigen::Matrix< double, 8, 1 >::Zero( ), 0 ).second.at( 0 );

        manualPartial.block( 0, j, 13, 1 ) =
                ( upPerturbedState.segment( 0, 13 ) - downPerturbedState.segment( 0, 13 ) ) / ( 2.0 * statePerturbation( j ) );

    }

    // Check element separately
    BOOST_CHECK_SMALL( std::fabs( manualPartial( 3, 1 + 10 ) -
                                  stateTransitionAndSensitivityMatrixAtEpoch( 3, 1 + 10 ) ), 1.0E-2 );
    manualPartial( 3, 1 + 10 ) = stateTransitionAndSensitivityMatrixAtEpoch( 3, 1 + 10 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualPartial.block( 0, 10, 6, 3 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 10, 6, 3 ) ), 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualPartial.block( 6, 10, 7, 3 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 6, 10, 7, 3 ) ), 1.0E-5 );

    // Numerically compute sensitivity matrix
    for( int j = 0; j < numberOfParametersToEstimate; j ++ )
    {
        Eigen::Matrix< double, 13, 1 > appliedStateDifference;

        Eigen::VectorXd upPerturbedState, downPerturbedState;
        perturbedState.setZero( );
        perturbedParameter.setZero( );
        perturbedParameter( j ) += parameterPerturbation( j );

        std::cout<<"Test "<<j<<" "<<perturbedParameter.transpose( )<<std::endl;

        upPerturbedState = executePhobosRotationSimulation< double, double >(
                    perturbedState, appliedStateDifference, perturbedParameter, 0 ).second.at( 0 );

        perturbedParameter.setZero( );
        perturbedParameter( j ) -= parameterPerturbation( j );
        downPerturbedState = executePhobosRotationSimulation< double, double >(
                    perturbedState, appliedStateDifference, perturbedParameter, 0 ).second.at( 0 );

        manualPartial.block( 0, j + 13, 13, 1 ) =
                ( upPerturbedState.segment( 0, 13 ) - downPerturbedState.segment( 0, 13 ) ) / ( 2.0 * parameterPerturbation( j ) );
    }
    std::cout<<manualPartial<<std::endl<<std::endl
            <<stateTransitionAndSensitivityMatrixAtEpoch<<std::endl<<std::endl<<
              ( manualPartial - stateTransitionAndSensitivityMatrixAtEpoch ).cwiseQuotient(
                  stateTransitionAndSensitivityMatrixAtEpoch )<<std::endl;

    // Check three values separately: could not find perturbations for which all partials are sufficiently within the linear regime

    BOOST_CHECK_SMALL( std::fabs( manualPartial( 4, 5 + 13 ) -
                                  stateTransitionAndSensitivityMatrixAtEpoch( 4, 5 + 13 ) ), 1.0E-4 );
    BOOST_CHECK_SMALL( std::fabs( manualPartial( 11, 4 + 13 ) -
                                  stateTransitionAndSensitivityMatrixAtEpoch( 11, 4 + 13 ) ), 1.0E-5 );
    BOOST_CHECK_SMALL( std::fabs( manualPartial( 12, 3 + 13 ) -
                                  stateTransitionAndSensitivityMatrixAtEpoch( 12, 3 + 13 ) ), 1.0E-2 );

    manualPartial( 4, 5 + 13 ) =  stateTransitionAndSensitivityMatrixAtEpoch( 4, 5 + 13 );
    manualPartial( 11, 4 + 13 ) =  stateTransitionAndSensitivityMatrixAtEpoch( 11, 4 + 13 );
    manualPartial( 12, 3 + 13 ) =  stateTransitionAndSensitivityMatrixAtEpoch( 12, 3 + 13 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( manualPartial.block( 0, 13, 13, 8 ) ), ( stateTransitionAndSensitivityMatrixAtEpoch.block( 0, 13, 13, 8 ) ), 1.0E-4 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

