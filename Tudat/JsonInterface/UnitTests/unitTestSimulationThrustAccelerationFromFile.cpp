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

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/jsonInterface.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_simulationThrustAccelerationFromFile )

BOOST_AUTO_TEST_CASE( test_json_simulationThrustAccelerationFromFile_main )
{

    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace numerical_integrators;
    using namespace interpolators;
    using namespace unit_conversions;
    using namespace propagators;
    using namespace json_interface;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             JSON SIMULATION              ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    JsonSimulationManager< > jsonSimulation( INPUT( "main" ) );
    jsonSimulation.updateSettings( );
    jsonSimulation.runPropagation( );
    std::map< double, Eigen::VectorXd > jsonResults =
            jsonSimulation.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > jsonResultsDependent =
            jsonSimulation.getDynamicsSimulator( )->getDependentVariableHistory( );



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

    // Define body settings for simulation.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = std::make_shared< BodySettings >( );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    double vehicleMass = 5000.0;
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define data to be used for thrust as a function of time.
    std::shared_ptr< FromFileDataMapSettings< Eigen::Vector3d > > thrustDataSettings =
            std::make_shared< FromFileDataMapSettings< Eigen::Vector3d > >( "thrustValues.txt" );

    // Define interpolator settings.
    std::shared_ptr< InterpolatorSettings >
            thrustInterpolatorSettings = std::make_shared< InterpolatorSettings >( linear_interpolator );

    // Create data interpolation settings
    std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolatorSettings =
            std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
                thrustDataSettings, thrustInterpolatorSettings );

    // Define specific impulse
    double constantSpecificImpulse = 3000.0;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back(
                std::make_shared< ThrustAccelerationSettings >(
                    thrustDataInterpolatorSettings, constantSpecificImpulse, lvlh_thrust_frame, "Earth" ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set initial conditions for the vehicle satellite that will be propagated in this simulation.
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );
    systemInitialState( 0 ) = 8.0E6;
    systemInitialState( 4 ) = 7.5E3;

    // Define propagation termination conditions (stop after 2 weeks).
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( 4.0E+5 );

    // Define settings for propagation of translational dynamics.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
              cowell );

    // Crete mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = createMassRateModel( "Vehicle", std::make_shared< FromThrustMassModelSettings >( 1 ),
                                                       bodyMap, accelerationModelMap );

    // Create settings for propagating the mass of the vehicle
    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Vehicle" }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                terminationSettings );

    // Create list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    lvlh_to_inertial_frame_rotation_dependent_variable, "Vehicle", "Earth" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );

    // Create propagation settings for mass and translational dynamics concurrently
    std::shared_ptr< PropagatorSettings< > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

    // Define integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 30.0 );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            //////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    const std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator =
            std::make_shared< SingleArcDynamicsSimulator< > >( bodyMap, integratorSettings, propagatorSettings );
    const std::map< double, Eigen::VectorXd > results = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
    const std::map< double, Eigen::VectorXd > resultsDependent = dynamicsSimulator->getDependentVariableHistory( );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        COMPARE RESULTS           ////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Check epochs, positions and velocities

    const std::vector< unsigned int > indices = { 0, 3, 6 };
    const std::vector< unsigned int > sizes = { 3, 3, 1 };
    const double tolerance = 1.0E-15;

    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonResults, results, indices, sizes, tolerance );


    // Check dependent variables

    const std::vector< unsigned int > indicesD = { 0, 3 };
    const std::vector< unsigned int > sizesD = { 3, 9 };

    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonResultsDependent, resultsDependent, indicesD, sizesD, tolerance );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////      CHECK CONSISTENCY OF from_json AND to_json FUNCTIONS          //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Convert jsonSimulation to JSON (using to_json functions) and use that to reset the simulation
    // (using from_json functions)
    jsonSimulation.resetJsonObject( jsonSimulation.getAsJson( ) );
    jsonSimulation.updateSettings( );

    // Get results
    jsonSimulation.runPropagation( );
    jsonResults = jsonSimulation.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
    jsonResultsDependent = jsonSimulation.getDynamicsSimulator( )->getDependentVariableHistory( );

    // Check epochs, positions and velocities
    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonResults, results, indices, sizes, tolerance );

    // Check dependent variables
    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonResultsDependent, resultsDependent, indicesD, sizesD, tolerance );

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
