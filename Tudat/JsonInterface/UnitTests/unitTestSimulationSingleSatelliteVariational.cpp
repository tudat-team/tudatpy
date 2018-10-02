/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include "Tudat/SimulationSetup/tudatEstimationHeader.h"
#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/jsonInterfaceVariational.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_simulationSingleSatelliteVariational )

BOOST_AUTO_TEST_CASE( test_json_simulationSingleSatelliteVariational_main )
{

    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace unit_conversions;
    using namespace json_interface;
    using namespace estimatable_parameters;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             JSON SIMULATION              ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    JsonVariationalEquationsSimulationManager< > jsonSimulation( INPUT( "main" ) );
    jsonSimulation.updateSettings( );
    jsonSimulation.runPropagation( );
    std::map< double, Eigen::VectorXd > jsonStates =
            jsonSimulation.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::MatrixXd > jsonStateTransition =
            jsonSimulation.getVariationalEquationsSolver( )->getNumericalVariationalEquationsSolution( )[ 0 ];
    std::map< double, Eigen::MatrixXd > jsonSensitivity =
            jsonSimulation.getVariationalEquationsSolver( )->getNumericalVariationalEquationsSolution( )[ 1 ];


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            MANUAL SIMULATION             ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationEndEpoch = 3600.0;


    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ) );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = 1.4888;
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 4.1137;
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.4084;
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = 2.4412;

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );


    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, fixedStepSize );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                    "Asterix", asterixInitialState, "Earth" ) );
    parameterNames.push_back(
                std::make_shared< EstimatableParameterSettings >(
                    "Earth", gravitational_parameter ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            //////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    const std::shared_ptr< SingleArcVariationalEquationsSolver< > > variationalEquationsSolver =
            std::make_shared< SingleArcVariationalEquationsSolver< > >(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, true, nullptr, false, true, false );

    const std::map< double, Eigen::VectorXd > states = variationalEquationsSolver->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    const std::map< double, Eigen::MatrixXd > stateTransition = variationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 0 ];
    const std::map< double, Eigen::MatrixXd > sensitivity = variationalEquationsSolver->getNumericalVariationalEquationsSolution( )[ 1 ];


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        COMPARE RESULTS           ////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const std::vector< unsigned int > indices = { 0, 3 };
    const std::vector< unsigned int > sizes = { 3, 3 };
    const double tolerance = 1.0E-15;

    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonStates, states, indices, sizes, tolerance );

    Eigen::MatrixXd finalStateTransition = jsonStateTransition.rbegin( )->second;

    std::vector< std::pair< unsigned int, unsigned int > > indicesStateTransition =
    { { 0, 0 }, { 0, 3 }, { 3, 0 }, { 3, 3 } };
    std::vector< std::pair< unsigned int, unsigned int > > sizesStateTransition =
    { { 3, 3 }, { 3, 3 }, { 3, 3 }, { 3, 3 } };

    std::vector< double > absoluteTolerancesStateTransition;
    absoluteTolerancesStateTransition.push_back(
                finalStateTransition.block( 0, 0, 3, 3 ).norm( ) * tolerance );
    absoluteTolerancesStateTransition.push_back(
                finalStateTransition.block( 0, 3, 3, 3 ).norm( ) * tolerance );
    absoluteTolerancesStateTransition.push_back(
                finalStateTransition.block( 3, 0, 3, 3 ).norm( ) * tolerance );
    absoluteTolerancesStateTransition.push_back(
                finalStateTransition.block( 3, 3, 3, 3 ).norm( ) * tolerance );

    BOOST_CHECK_CLOSE_INTEGRATION_MATRIX_RESULTS(
                stateTransition, jsonStateTransition, indicesStateTransition,
                sizesStateTransition, absoluteTolerancesStateTransition );



    Eigen::MatrixXd jsonFinalSensitivity = jsonSensitivity.rbegin( )->second;

    for( int i = 0; i < jsonFinalSensitivity.cols( ); i++ )
    {
        Eigen::MatrixXd jsonFinalSensitivity = jsonSensitivity.rbegin( )->second;

        std::vector< std::pair< unsigned int, unsigned int > > indicesSensitivity;
        std::vector< std::pair< unsigned int, unsigned int > > sizesSensitivity;
        std::vector< double > absoluteTolerancesSensitivity;

        indicesSensitivity.push_back( std::make_pair( 0, i ) );
        indicesSensitivity.push_back( std::make_pair( 3, i ) );

        sizesSensitivity.push_back( std::make_pair( 3, 1 ) );
        sizesSensitivity.push_back( std::make_pair( 3, 1 ) );

        absoluteTolerancesSensitivity.push_back(
                    jsonFinalSensitivity.block( 0, i, 3, 1 ).norm( ) * std::numeric_limits< double >::epsilon( ) );
        absoluteTolerancesSensitivity.push_back(
                    jsonFinalSensitivity.block( 3, i, 3, 1 ).norm( ) * std::numeric_limits< double >::epsilon( ) );


        BOOST_CHECK_CLOSE_INTEGRATION_MATRIX_RESULTS(
                    sensitivity, jsonSensitivity, indicesSensitivity,
                    sizesSensitivity, absoluteTolerancesSensitivity );
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////      CHECK CONSISTENCY OF from_json AND to_json FUNCTIONS          //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Convert jsonSimulation to JSON (using to_json functions) and use that to reset the simulation
//    // (using from_json functions)
//    jsonSimulation.resetJsonObject( jsonSimulation.getAsJson( ) );
//    jsonSimulation.updateSettings( );

//    // Get results
//    jsonSimulation.runPropagation( );
//    jsonStates =
//            jsonSimulation.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
//    jsonStateTransition =
//            jsonSimulation.getVariationalEquationsSolver( )->getNumericalVariationalEquationsSolution( )[ 0 ];
//    jsonSensitivity =
//            jsonSimulation.getVariationalEquationsSolver( )->getNumericalVariationalEquationsSolution( )[ 1 ];

//    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonStates, states, indices, sizes, tolerance );

//    BOOST_CHECK_CLOSE_INTEGRATION_MATRIX_RESULTS(
//                stateTransition, jsonStateTransition, indicesStateTransition,
//                sizesStateTransition, absoluteTolerancesStateTransition );

//    BOOST_CHECK_CLOSE_INTEGRATION_MATRIX_RESULTS(
//                sensitivity, jsonSensitivity, indicesSensitivity,
//                sizesSensitivity, absoluteTolerancesSensitivity );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
