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

BOOST_AUTO_TEST_SUITE( test_json_simulationGalileoConstellation )

BOOST_AUTO_TEST_CASE( test_json_simulationGalileoConstellation_main )
{

    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace basic_astrodynamics;
    using namespace basic_mathematics;
    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using namespace input_output;
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



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            MANUAL SIMULATION             ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLES      ////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3600.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 30.0;

    // Set number of satellites in constellation.
    const unsigned int numberOfSatellites = 30;

    // Set number of planes in constellation.
    const unsigned int numberOfPlanes = 3;

    // Set number of satellites per plane in constellation.
    const unsigned int numberOfSatellitesPerPlane = numberOfSatellites / numberOfPlanes;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define environment settings
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings = getDefaultBodySettings( { "Earth" } );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = nullptr;
    bodySettings[ "Earth" ]->shapeModelSettings = nullptr;

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    // Set accelerations for each satellite.
    std::string currentSatelliteName;
    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        currentSatelliteName =  "Satellite" + std::to_string( i );
        bodyMap[ currentSatelliteName ] = std::make_shared< simulation_setup::Body >( );
    }

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE CONSTELLATION INITIAL STATES      /////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set orbital parameters of Galileo constellation.
    const double semiMajorAxis = 23222.0e3 + 6378.1e3;                                // [km]
    const double eccentricity = 0.0;                                                  // [-]
    const double inclination = convertDegreesToRadians( 56.0 );                       // [rad]
    const double argumentOfPeriapsis = 0.0;                                           // [rad]
    const double longitudeOfAscendingNodeSpacing
            = 2.0 * mathematical_constants::PI / numberOfPlanes;                          // [rad]
    const double trueAnomalySpacing
            = 2.0 * mathematical_constants::PI / numberOfSatellitesPerPlane;              // [rad]

    // Set initial conditions for Galileo satellites that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to
    // Cartesian elements. They are stored in a matrix.

    // Declare size of state.
    const unsigned int sizeOfState = 6;

    // Set Keplerian elements for Galileo satellites.
    Eigen::MatrixXd initialConditionsInKeplerianElements( sizeOfState, numberOfSatellites );

    // Set semiMajorAxis.
    initialConditionsInKeplerianElements.row( 0 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, semiMajorAxis );

    // Set eccentricity.
    initialConditionsInKeplerianElements.row( 1 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, eccentricity );

    // Set inclination.
    initialConditionsInKeplerianElements.row( 2 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, inclination );

    // Set argument of periapsis.
    initialConditionsInKeplerianElements.row( 3 )
            = Eigen::MatrixXd::Constant( 1, numberOfSatellites, argumentOfPeriapsis );

    // Set longitude of ascending node.
    for ( unsigned int i = 0; i < numberOfPlanes; i++ )

    {
        initialConditionsInKeplerianElements.block( 4, i * numberOfSatellitesPerPlane,
                                                    1, numberOfSatellitesPerPlane )
                = Eigen::MatrixXd::Constant( 1, numberOfSatellitesPerPlane,
                                             i * longitudeOfAscendingNodeSpacing );
    }

    // Set true anomaly.
    Eigen::RowVectorXd trueAnomalySpacingIntegers( numberOfSatellitesPerPlane );

    for ( unsigned int i = 0; i < numberOfSatellitesPerPlane; i++ )
    {
        trueAnomalySpacingIntegers( i ) =  i * 1.0;
    }

    for ( unsigned int i = 0; i < numberOfPlanes; i++ )
    {
        initialConditionsInKeplerianElements.block( 5, i * numberOfSatellitesPerPlane,
                                                    1, numberOfSatellitesPerPlane )
                = Eigen::MatrixXd::Constant( 1, numberOfSatellitesPerPlane,
                                             trueAnomalySpacing ).array( )
                * trueAnomalySpacingIntegers.array( );
    }

    // Convert initial conditions to Cartesian elements.
    Eigen::MatrixXd initialConditions( sizeOfState, numberOfSatellites );

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        Eigen::Vector6d initKepl = initialConditionsInKeplerianElements.col( i ).cast< double >();
        initialConditions.col( i ) = convertKeplerianToCartesianElements(
                    initKepl, static_cast< double >(earthGravitationalParameter) );
    }

    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 * numberOfSatellites );
    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        systemInitialState.segment( i * 6, 6 ) = initialConditions.col( i );
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            /////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Set accelerations for each satellite.
    for ( unsigned int i = 0; i < numberOfSatellites; i++ )
    {
        currentSatelliteName =  "Satellite" + std::to_string( i );

        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCurrentSatellite;
        accelerationsOfCurrentSatellite[ "Earth" ].push_back(
                    std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
        accelerationMap[ currentSatelliteName ] = accelerationsOfCurrentSatellite;

        bodiesToPropagate.push_back( currentSatelliteName );
        centralBodies.push_back( "Earth" );
    }


    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            //////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    const std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator =
            std::make_shared< SingleArcDynamicsSimulator< > >( bodyMap, integratorSettings, propagatorSettings );
    const std::map< double, Eigen::VectorXd > results = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        COMPARE RESULTS           ////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< unsigned int > indices;
    std::vector< unsigned int > sizes;
    for ( unsigned int i = 0; i < numberOfSatellites; ++i )
    {
        indices.push_back( 3 * i );
        sizes.push_back( 3 );
    }
    const double tolerance = 1.0E-15;

    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonResults, results, indices, sizes, tolerance );



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

    BOOST_CHECK_CLOSE_INTEGRATION_RESULTS( jsonResults, results, indices, sizes, tolerance );

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
