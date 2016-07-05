/*    Copyright (c) 2010-2016, Delft University of Technology
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
#include <omp.h>

#include <boost/make_shared.hpp>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_encke_propagator )

// Test Encke propagator for point mass central body.
BOOST_AUTO_TEST_CASE( testEnckePopagatorForPointMassCentralBodies )
{
    // Test simulation for different central body cases
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
        std::string kernelsPath = input_output::getSpiceKernelPath( );
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
        spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

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
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
        accelerationsOfEarth[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfEarth[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfEarth[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Earth" ] = accelerationsOfEarth;

        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
        accelerationsOfMars[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMars[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMars[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Mars" ] = accelerationsOfMars;

        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
        accelerationsOfMoon[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfMoon[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
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
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4,
                  initialEphemerisTime, 250.0 );

        // Create propagation settings (Cowell)
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime );

        // Propagate orbit with Cowell method
        SingleArcDynamicsSimulator< double > dynamicsSimulator2(
                    bodyMap, integratorSettings, propagatorSettings, true );

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
        propagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime, encke );

        // Propagate orbit with Encke method
        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true );

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
        spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
        spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
        spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

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
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        }
        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Create spacecraft object.
        bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );
        bodyMap[ "Vehicle" ]->setEphemeris( boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                boost::shared_ptr< interpolators::OneDimensionalInterpolator
                                                < double, Vector6d  > >( ), "Earth", "J2000" ) );
        boost::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings =
                boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", 4.0, 1.2, boost::assign::list_of( "Earth" )( "Moon" ) );
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
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

        // Use only central gravity for Earth
        if( simulationCase < 2 )
        {
            accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
        }
        // Use spherical harmonics for Earth
        else
        {
            accelerationsOfVehicle[ "Earth" ].push_back(
                        boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

        }

        // Use perturbations other than Earth gravity
        if( simulationCase % 2 == 0 )
        {
            accelerationsOfVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Venus" ].push_back( boost::make_shared< AccelerationSettings >(
                                                             basic_astrodynamics::central_gravity ) );
            accelerationsOfVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::cannon_ball_radiation_pressure ) );
        }
        accelerationMap[  "Vehicle" ] = accelerationsOfVehicle;
        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        // Set Keplerian elements for Vehicle.
        Vector6d vehicleInitialStateInKeplerianElements;
        vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 8000.0E3;
        vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
        vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
        vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        const Vector6d vehicleInitialState = convertKeplerianToCartesianElements(
                    vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

        // Define propagator settings (Cowell)
        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, simulationEndEpoch );

        // Define integrator settings.
        const double fixedStepSize = 5.0;
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, fixedStepSize );

        // Propagate orbit with Cowell method
        SingleArcDynamicsSimulator< double > dynamicsSimulator2(
                    bodyMap, integratorSettings, propagatorSettings, true );

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
                    bodyMap[ "Vehicle" ]->getEphemeris( )->getCartesianStateFromEphemeris( currentTestTime );

            currentTestTime += testTimeStep;
        }

        // Create propagation settings (Encke)
        propagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, simulationEndEpoch, encke );

        // Propagate orbit with Encke method
        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true );

        // Get resutls of Encke integration at given times.
        currentTestTime = initialTestTime;
        std::map< double, Eigen::Matrix< double, 6, 1 > > enckeIntegrationResults;
        while( currentTestTime < finalTestTime )
        {
            enckeIntegrationResults[ currentTestTime ].segment( 0, 6 ) =
                    bodyMap[ "Vehicle" ]->getEphemeris( )->getCartesianStateFromEphemeris( currentTestTime );
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
BOOST_AUTO_TEST_SUITE_END( )


}

}


