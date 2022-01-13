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

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_arcwise_environment )

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


//! Unit test to check if tidal time lag parameters are estimated correctly
BOOST_AUTO_TEST_CASE( test_ArcwiseEnvironmentParameters )
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

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.01, 0.1 ).finished( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
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
    bodies.at( "Vehicle" )->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Vehicle", bodies ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::point_mass_gravity ) );
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
    std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable, "Vehicle", "Earth" ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    radiation_pressure_coefficient_dependent_variable, "Vehicle", "Sun" ) );
    dependentVariableSaveSettings = std::make_shared< DependentVariableSaveSettings >( dependentVariables );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, double( finalEphemerisTime ),
                cowell, dependentVariableSaveSettings);

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                double( initialEphemerisTime ), 90.0,
                RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
                90.0, 90.0, 1.0, 1.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE OBSERVATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    LinkEnds linkEnds;
    linkEnds[ observed_body ] = std::make_pair( "Vehicle", "" );
    std::vector< std::shared_ptr< ObservationModelSettings > >  observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnds ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );

    std::vector< double > arcStartTimeList =
    { initialEphemerisTime, initialEphemerisTime + 0.33 * 86400.0, initialEphemerisTime + 0.66 * 86400.0 };
    parameterNames.push_back(
                std::make_shared< ArcWiseDragCoefficientEstimatableParameterSettings >( "Vehicle", arcStartTimeList ) );
    parameterNames.push_back(
                std::make_shared< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings >( "Vehicle", arcStartTimeList ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );

    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    truthParameters( 7 ) += 0.1;
    truthParameters( 8 ) += 0.2;

    truthParameters( 10 ) += 0.1;
    truthParameters( 11 ) += 0.2;

    parametersToEstimate->resetParameterValues( truthParameters );
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          INITIALIZE ORBIT DETERMINATION OBJECT     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );

    std::map< double, Eigen::VectorXd > dependentVariableData =
            orbitDeterminationManager.getVariationalEquationsSolver( )->getDynamicsSimulatorBase(
                )->getDependentVariableNumericalSolutionBase( )[ 0 ];

    // Test whether arc-wise coefficients are correctly used.
    double testDragCoefficient = 0.0;
    double testRadiationPressureCoefficient = 0.0;
    for( auto variableIterator : dependentVariableData )
    {
        double currentTime = variableIterator.first;

        if( currentTime < arcStartTimeList.at( 1 ) )
        {
            testDragCoefficient = truthParameters( 6 );
            testRadiationPressureCoefficient = truthParameters( 9 );
        }
        else if( currentTime < arcStartTimeList.at( 2 ) )
        {
            testDragCoefficient = truthParameters( 7 );
            testRadiationPressureCoefficient = truthParameters( 10 );
        }
        else
        {
            testDragCoefficient = truthParameters( 8 );
            testRadiationPressureCoefficient = truthParameters( 11 );
        }

        BOOST_CHECK_EQUAL( testDragCoefficient, variableIterator.second( 0 ) );
        BOOST_CHECK_EQUAL( aerodynamicCoefficient * -0.01, variableIterator.second( 1 ) );
        BOOST_CHECK_EQUAL( aerodynamicCoefficient * 0.1, variableIterator.second( 2 ) );
        BOOST_CHECK_EQUAL( testRadiationPressureCoefficient, variableIterator.second( 3 ) );

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          SIMULATE OBSERVATIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define time of first observation
    double observationTimeStart = initialEphemerisTime + 1000.0;

    // Define time between two observations
    double  observationInterval = 30.0;

    // Simulate observations for 3 days
    std::vector< double > baseTimeList;
    double currentTime = initialEphemerisTime + 3600.0;
    while( currentTime < finalEphemerisTime - 3600.0 )
    {
        baseTimeList.push_back( currentTime );
        currentTime += observationInterval;
    }

    // Create measureement simulation input
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< > >(
                    position_observable, linkEnds, baseTimeList, observed_body ) );

    // Simulate observations
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 10.0 );
    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.0E-2 );
    parameterPerturbation.segment( 6, 6 ) = Eigen::VectorXd::Constant( 6, 1.0 );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate = truthParameters;
    initialParameterEstimate += parameterPerturbation;


    // Define estimation input
    std::shared_ptr< PodInput< double, double > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, initialParameterEstimate.rows( ),
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
                initialParameterEstimate - truthParameters );
    podInput->defineEstimationSettings( true, true, false, true );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 4 ) );
    Eigen::VectorXd parameterEstimate = podOutput->parameterEstimate_ - truthParameters;

    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( parameterEstimate( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( parameterEstimate( i + 3 ) ), 1.0E-4 );
        BOOST_CHECK_SMALL( std::fabs( parameterEstimate( i + 6 ) ), 1.0E-5 );
        BOOST_CHECK_SMALL( std::fabs( parameterEstimate( i + 9 ) ), 1.0E-4 );


    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
