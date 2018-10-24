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

#include <limits>
#include <string>
#include <thread>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EstimationSetup/createNumericalSimulator.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_encke_propagator )

// Test Encke propagator for point mass central body.
BOOST_AUTO_TEST_CASE( testEnckePopagatorForPointMassCentralBodies )
{
    // Test simulation for diffe  rent central body cases
    for( unsigned int simulationCase = 0; simulationCase < 2; simulationCase++ )
    {
        //Using declarations.
        using namespace tudat::interpolators;
        using namespace tudat::numerical_integrators;
        using namespace tudat::spice_interface;
        using namespace tudat::simulation_setup;
        using namespace tudat::basic_astrodynamics;
        using namespace tudat::orbital_element_conversions;
        using namespace tudat::propagators;


        //Load spice kernels.
        spice_interface::loadStandardSpiceKernels( );

        // Define bodies in simulation.
        unsigned int totalNumberOfBodies = 7;
        std::vector< std::string > bodyNames;
        bodyNames.resize( totalNumberOfBodies );
        bodyNames[ 0 ] = "Earth";
        bodyNames[ 1 ] = "Mars";
        bodyNames[ 2 ] = "Sun";
        bodyNames[ 3 ] = "Venus";
        bodyNames[ 4 ] = "Moon";
        bodyNames[ 5 ] = "Mercury";
        bodyNames[ 6 ] = "Jupiter";

        double initialEphemerisTime = 1.0E7;
        double finalEphemerisTime = 2.0E7;
        double maximumTimeStep = 3600.0;
        double buffer = 5.0 * maximumTimeStep;

        // Create bodies needed in simulation
        NamedBodyMap bodyMap = createBodies(
                    getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer ) );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
        accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfEarth[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Earth" ] = accelerationsOfEarth;

        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
        accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMars[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Mars" ] = accelerationsOfMars;

        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
        accelerationsOfMoon[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Moon" ] = accelerationsOfMoon;

        // Propagate Earth, Mars and Moon
        std::vector< std::string > bodiesToPropagate;
        bodiesToPropagate.push_back( "Earth" );
        bodiesToPropagate.push_back( "Mars" );
        bodiesToPropagate.push_back( "Moon" );

        unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );

        // Define central bodies: all Sun for simulationCase = 0, Earth and Mars: Sun, Moon: Earth for simulationCase = 1
        std::vector< std::string > centralBodies;
        std::map< std::string, std::string > centralBodyMap;
        centralBodies.resize( numberOfNumericalBodies );
        for( int i = 0; i < 3; i++ )
        {
            if( i == 2 && simulationCase == 1 )
            {
                centralBodies[ i ] = "Earth";
            }
            else
            {
                centralBodies[ i ] = "Sun";
            }
            centralBodyMap[ bodiesToPropagate[ i ] ] = centralBodies[ i ];
        }


        // Get initial states for bodies.
        Eigen::VectorXd systemInitialState = Eigen::VectorXd( bodiesToPropagate.size( ) * 6 );
        for( unsigned int i = 0; i < numberOfNumericalBodies ; i++ )
        {
            systemInitialState.segment( i * 6 , 6 ) =
                    bodyMap[ bodiesToPropagate[ i ] ]->getStateInBaseFrameFromEphemeris( initialEphemerisTime ) -
                    bodyMap[ centralBodies[ i ] ]->getStateInBaseFrameFromEphemeris( initialEphemerisTime );
        }

        // Create acceleratiuon models.
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, centralBodyMap );

        // Create integrator settings.
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4,
                  initialEphemerisTime, 250.0 );

        // Create propagation settings (Cowell)
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime );

        // Propagate orbit with Cowell method
        SingleArcDynamicsSimulator< double > dynamicsSimulator2(
                    bodyMap, integratorSettings, propagatorSettings, true, false, true );

        // Define ephemeris interrogation settings.
        double initialTestTime = initialEphemerisTime + 10.0 * maximumTimeStep;
        double finalTestTime = finalEphemerisTime - 10.0 * maximumTimeStep;
        double testTimeStep = 1.0E4;

        // Get resutls of Cowell integration at given times.
        double currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 18, 1 > > cowellIntegrationResults;
        while( currentTestTime < finalTestTime )
        {
            cowellIntegrationResults[ currentTestTime ].segment( 0, 6 ) =
                    bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            cowellIntegrationResults[ currentTestTime ].segment( 6, 6 ) =
                    bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            cowellIntegrationResults[ currentTestTime ].segment( 12, 6 ) =
                    bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( currentTestTime );

            currentTestTime += testTimeStep;
        }

        // Create propagation settings (Encke)
        propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime, encke );

        // Propagate orbit with Encke method
        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, true );

        // Get resutls of Encke integration at given times.
        currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 18, 1 > > enckeIntegrationResults;
        while( currentTestTime < finalTestTime )
        {
            enckeIntegrationResults[ currentTestTime ].segment( 0, 6 ) =
                    bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            enckeIntegrationResults[ currentTestTime ].segment( 6, 6 ) =
                    bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            enckeIntegrationResults[ currentTestTime ].segment( 12, 6 ) =
                    bodyMap[ "Moon" ]->getStateInBaseFrameFromEphemeris( currentTestTime );
            currentTestTime += testTimeStep;
        }

        // Compare results of Cowell and Encke propagations
        std::map< double, Eigen::Matrix< double, 18, 1 > >::iterator enckeIterator = enckeIntegrationResults.begin( );
        std::map< double, Eigen::Matrix< double, 18, 1 > >::iterator cowellIterator = cowellIntegrationResults.begin( );
        for( unsigned int i = 0; i < enckeIntegrationResults.size( ); i++ )
        {
            for( int j= 0; j< 3; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second ).segment( j, 1 )( 0 ), 0.01 );
            }

            for( int j = 6; j < 9; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second ).segment( j, 1 )( 0 ), 0.01 );
            }

            for( int j = 12; j < 15; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second ).segment( j, 1 )( 0 ), 0.075 );
            }

            for( int j = 3; j < 6; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second ).segment( j, 1 )( 0 ), 1.0E-8 );
            }

            for( int j = 9; j < 12; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second ).segment( j, 1 )( 0 ), 1.0E-8 );

            }

            for( int j = 15; j < 18; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second ).segment( j, 1 )( 0 ), 1.0E-6 );

            }
            enckeIterator++;
            cowellIterator++;
        }
    }
}

// Test Encke propagator for point mass, and spherical harmonics central body.
BOOST_AUTO_TEST_CASE( testEnckePopagatorForSphericalHarmonicCentralBodies )
{
    for( unsigned int simulationCase = 0; simulationCase < 4; simulationCase++ )
    {
        using namespace tudat;
        using namespace simulation_setup;
        using namespace propagators;
        using namespace numerical_integrators;
        using namespace orbital_element_conversions;
        using namespace basic_mathematics;
        using namespace gravitation;
        using namespace numerical_integrators;

        // Load Spice kernels.
        spice_interface::loadStandardSpiceKernels( );

        // Set simulation time settings.
        const double simulationStartEpoch = 0.0;
        const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );
        bodiesToCreate.push_back( "Mars" );
        bodiesToCreate.push_back( "Venus" );

        // Create body objects.
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        }
        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Create spacecraft object.
        bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
        bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );
        bodyMap[ "Vehicle" ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                < double, Eigen::Vector6d  > >( ), "Earth", "J2000" ) );
        std::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings =
                std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", 4.0, 1.2, std::vector< std::string >{ "Earth", "Moon" } );
        bodyMap[ "Vehicle" ]->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        vehicleRadiationPressureSettings, "Vehicle", bodyMap ) );


        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

        // Use only central gravity for Earth
        if( simulationCase < 2 )
        {
            accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
        }
        // Use spherical harmonics for Earth
        else
        {
            accelerationsOfVehicle[ "Earth" ].push_back(
                        std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

        }

        // Use perturbations other than Earth gravity
        if( simulationCase % 2 == 0 )
        {
            accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::cannon_ball_radiation_pressure ) );
        }
        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Set Keplerian elements for Vehicle.
        Eigen::Vector6d vehicleInitialStateInKeplerianElements;
        vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 8000.0E3;
        vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
        vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
        vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        const Eigen::Vector6d vehicleInitialState = convertKeplerianToCartesianElements(
                    vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

        // Define propagator settings (Cowell)
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, simulationEndEpoch );

        // Define integrator settings.
        const double fixedStepSize = 5.0;
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, fixedStepSize );

        // Propagate orbit with Cowell method
        SingleArcDynamicsSimulator< double > dynamicsSimulator2(
                    bodyMap, integratorSettings, propagatorSettings, true, false, true );

        // Define ephemeris interrogation settings.
        double initialTestTime = simulationStartEpoch + 10.0 * fixedStepSize;
        double finalTestTime = simulationEndEpoch - 10.0 * fixedStepSize;
        double testTimeStep = 1.0E4;

        // Get resutls of Cowell integration at given times.
        double currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 6, 1 > > cowellIntegrationResults;
        while( currentTestTime < finalTestTime )
        {
            cowellIntegrationResults[ currentTestTime ].segment( 0, 6 ) =
                    bodyMap[ "Vehicle" ]->getEphemeris( )->getCartesianState( currentTestTime );

            currentTestTime += testTimeStep;
        }

        // Create propagation settings (Encke)
        propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, simulationEndEpoch, encke );

        // Propagate orbit with Encke method
        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, true );

        // Get resutls of Encke integration at given times.
        currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 6, 1 > > enckeIntegrationResults;
        while( currentTestTime < finalTestTime )
        {
            enckeIntegrationResults[ currentTestTime ].segment( 0, 6 ) =
                    bodyMap[ "Vehicle" ]->getEphemeris( )->getCartesianState( currentTestTime );
            currentTestTime += testTimeStep;
        }

        // Compare results of Cowell and Encke propagations
        std::map< double, Eigen::Matrix< double, 6, 1 > >::iterator enckeIterator = enckeIntegrationResults.begin( );
        std::map< double, Eigen::Matrix< double, 6, 1 > >::iterator cowellIterator = cowellIntegrationResults.begin( );
        for( unsigned int i = 0; i < enckeIntegrationResults.size( ); i++ )
        {
            for( int j= 0; j< 3; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second )( j ), 0.02 );
            }

            for( int j = 3; j < 6; j++ )
            {
                BOOST_CHECK_SMALL( ( enckeIterator->second - cowellIterator->second )( j ), 1.0E-5 );

            }
            enckeIterator++;
            cowellIterator++;
        }
    }
}

//! test if Encke propagator works properly for high eccentricities. Test if the propagation continmues unimpeded for
//! large eccentricities (was found to fail in rpevious version due to no convergence of root finder
BOOST_AUTO_TEST_CASE( testEnckePopagatorForHighEccentricities )
{
    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace unit_conversions;

    std::vector< double > testEccentricities;
    testEccentricities.push_back( 0.8 );
    testEccentricities.push_back( 0.9 );
    testEccentricities.push_back( 0.95 );
    testEccentricities.push_back( 0.99 );

    for( unsigned testCase = 0; testCase < testEccentricities.size( ); testCase++ )
    {
        // Set simulation end epoch.
        const double simulationEndEpoch = 2.0 * tudat::physical_constants::JULIAN_DAY;

        // Set numerical integration fixed step size.
        const double fixedStepSize = 15.0;

        // Define body settings for simulation.
        std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
        bodySettings[ "Earth" ] = std::make_shared< BodySettings >( );
        bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ), "SSB", "J2000" );
        bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

        // Create Earth object
        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Create spacecraft object.
        bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );


        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );


        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::central_gravity ) );
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Set Keplerian elements for Asterix.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = testEccentricities.at( testCase );
        asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements,
                    earthGravitationalParameter );

        TranslationalPropagatorType propagatorType = encke;
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, propagatorType);

        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( 0.0, fixedStepSize,
                  RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0E-4, 3600.0, 1.0E-14, 1.0E-14 );


        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        // Check if orbit is properly propagated
        Eigen::VectorXd analyticalSolution;
        double finalTime = integrationResult.rbegin( )->first;

        // Check if final time is correctly reached
        BOOST_CHECK_EQUAL( simulationEndEpoch <=  finalTime, true );

        // Loosen velocity tolerance for high eccentricity
        double toleranceMultiplier = ( testCase > 2 ) ? 100.0 : 1.0;

        // Compare propagated orbit against numerical result.
        for( std::map< double, Eigen::VectorXd >::const_iterator resultIterator = integrationResult.begin( );
             resultIterator != integrationResult.end( ); resultIterator++ )
        {
            analyticalSolution = convertKeplerianToCartesianElements( propagateKeplerOrbit(
                        asterixInitialStateInKeplerianElements, resultIterator->first, earthGravitationalParameter ),
                                         earthGravitationalParameter );
            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( analyticalSolution( i ) - resultIterator->second( i ) ), 1.0E-2 );
                BOOST_CHECK_SMALL( std::fabs( analyticalSolution( i  + 3 ) - resultIterator->second( i + 3  ) ),
                                   toleranceMultiplier * 5.0E-6 );

            }
        }
    }

}
BOOST_AUTO_TEST_SUITE_END( )


}

}


