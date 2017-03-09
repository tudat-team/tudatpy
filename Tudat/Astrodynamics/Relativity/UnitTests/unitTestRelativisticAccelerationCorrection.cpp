/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#define BOOST_TEST_MAIN

#include <boost/assign/list_of.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

namespace tudat
{

namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::unit_conversions;

double computeLenseThirringNodePrecession(
        const double gravitationalParameter,
        const double angularMomentum,
        const double semiMajorAxis,
        const double eccentricity,
        const double inclination )
{
    return - 6.0 * gravitationalParameter * angularMomentum * std::cos(inclination  ) /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * semiMajorAxis * semiMajorAxis * semiMajorAxis *
              std::pow( 1.0 - eccentricity * eccentricity, 1.5 ) );
}

double computeLenseThirringPericenterPrecession(
        const double gravitationalParameter,
        const double angularMomentum,
        const double semiMajorAxis,
        const double eccentricity )
{
    return 2.0 * gravitationalParameter * angularMomentum /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * semiMajorAxis * semiMajorAxis * semiMajorAxis *
              std::pow( 1.0 - eccentricity * eccentricity, 1.5 ) );
}

BOOST_AUTO_TEST_SUITE( test_relativistic_acceleration_corrections )

BOOST_AUTO_TEST_CASE( testLenseThirring )
{
    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation end epoch.
    const double simulationEndEpoch = 50.0 * tudat::physical_constants::JULIAN_DAY;


    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ) );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                     false, true, false ) );
    accelerationMap[  "Asterix" ] = accelerationsOfAsterix;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 5000.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.2;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 65.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, encke );


    // Create numerical integrator.
    double simulationStartEpoch = 0.0;
    const double fixedStepSize = 10.0;
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( rungeKuttaVariableStepSize, 0.0, 10.0,
              RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0E-3, 1.0E3, 1.0E-14, 1.0E-14 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > keplerianIntegrationResult;

    // Compute map of Kepler elements
    Eigen::Vector6d currentCartesianState;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
        currentCartesianState = stateIterator->second;
        keplerianIntegrationResult[ stateIterator->first ] =
                convertCartesianToKeplerianElements(
                    currentCartesianState, earthGravitationalParameter );
        //std::cout<<std::setprecision( 16 )<<stateIterator->first<<" "<<keplerianIntegrationResult[ stateIterator->first ].transpose( )<<std::endl;
    }

    std::cout<< computeLenseThirringPericenterPrecession(
                    earthGravitationalParameter, 1.0E9,  asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                    asterixInitialStateInKeplerianElements( eccentricityIndex ) )<<
                " "<< computeLenseThirringNodePrecession(
                    earthGravitationalParameter, 1.0E9,  asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                    asterixInitialStateInKeplerianElements( eccentricityIndex ),
                    asterixInitialStateInKeplerianElements( inclinationIndex ) )<<std::endl;


    //input_output::writeDataMapToTextFile( keplerianIntegrationResult, "ltKepler2.dat" );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
