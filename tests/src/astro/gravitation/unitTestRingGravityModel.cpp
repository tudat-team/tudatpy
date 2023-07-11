/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace propagators;
using namespace numerical_integrators;
using namespace orbital_element_conversions;
using namespace basic_mathematics;
using namespace gravitation;
using namespace numerical_integrators;
using namespace basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_ring_gravity_model )

BOOST_AUTO_TEST_CASE( testRingVersusPointMassesGravityModelSinglePoints )
{
    double tolerance = 1e-14;

    // Define ring gravity parameters
//    const double gravitationalParameter = 2.39e21 * physical_constants::GRAVITATIONAL_CONSTANT;
//    const double ringRadius = 2.7 * physical_constants::ASTRONOMICAL_UNIT;
    const double gravitationalParameter = 2.0;
    const double ringRadius = 1.5;

    // Number of point masses in discrete ring model
    int numPointMasses = 100;

    Eigen::Vector3d testPosition;

    for ( unsigned int i = 0; i < 7; ++i )
    {
        if ( i == 0 )
        {
            testPosition = ( Eigen::Vector3d( ) << 3.6, 0.0, 0.0 ).finished( );
        }
        else if ( i == 1 )
        {
            testPosition = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );
        }
        else if ( i == 2 )
        {
            testPosition = ( Eigen::Vector3d( ) << 0.0, 0.0, 1.0 ).finished( );
        }
        else if ( i == 3 )
        {
            testPosition = ( Eigen::Vector3d( ) << 0.2, 1.0, 0.0 ).finished( );
        }
        else if ( i == 4 )
        {
            testPosition = ( Eigen::Vector3d( ) << -4.0, 2.6, 0.0 ).finished( );
        }
        else if ( i == 5 )
        {
            testPosition = ( Eigen::Vector3d( ) << 0.6, 1.0, 0.4 ).finished( );
        }
        else if ( i == 6 )
        {
            testPosition = ( Eigen::Vector3d( ) << 0.7, -0.4, -0.2 ).finished( );
        }

        Eigen::Vector3d ringAcceleration = Eigen::Vector3d::Zero( ), pointMassAcceleration = Eigen::Vector3d::Zero( );
        double ringPotential = 0.0, pointMassPotential = 0.0;

        for ( unsigned int gravityModelId : { 0, 1 } )
        {
            std::function< void ( Eigen::Vector3d& ) > positionFunction = [ = ] ( Eigen::Vector3d& positionInput ) {
                positionInput = testPosition; };

            if ( gravityModelId == 0 )
            {
                RingGravitationalAccelerationModel ringGravityModel = RingGravitationalAccelerationModel(
                        positionFunction, gravitationalParameter, ringRadius, false );
                ringGravityModel.resetUpdatePotential( true );
                ringGravityModel.updateMembers( TUDAT_NAN );

                ringAcceleration = ringGravityModel.getAcceleration( );
                ringPotential = ringGravityModel.getCurrentPotential( );
            }
            else
            {
                for ( int j = 0; j < numPointMasses; ++j )
                {
                    Eigen::Vector3d pointMassPosition = Eigen::Vector3d::Zero( );
                    pointMassPosition( 0 ) = ringRadius * std::cos( j * 2.0 * mathematical_constants::PI / numPointMasses );
                    pointMassPosition( 1 ) = ringRadius * std::sin( j * 2.0 * mathematical_constants::PI / numPointMasses );

                    std::function< void ( Eigen::Vector3d& ) > pointMassPositionFunction = [ = ] (
                            Eigen::Vector3d& positionInput ) { positionInput = pointMassPosition; };

                    CentralGravitationalAccelerationModel3d pointMassGravityModel = CentralGravitationalAccelerationModel3d(
                            positionFunction, gravitationalParameter / numPointMasses, pointMassPositionFunction );
                    pointMassGravityModel.resetUpdatePotential( true );
                    pointMassGravityModel.updateMembers( TUDAT_NAN );

                    pointMassPotential += pointMassGravityModel.getCurrentPotential( );
                    pointMassAcceleration += pointMassGravityModel.getAcceleration( );
                }
            }
        }

        // For cases in which acceleration is 0, add 1 (otherwise relative error test will fail)
        for ( unsigned int i = 0; i < ringAcceleration.size( ); ++i )
        {
            if ( std::abs( ringAcceleration( i ) ) < 1e-16 )
            {
                ringAcceleration( i ) += 1;
                pointMassAcceleration( i ) += 1;
            }
        }

        BOOST_CHECK_CLOSE_FRACTION( pointMassPotential, ringPotential, tolerance );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( pointMassAcceleration, ringAcceleration, tolerance );
    }

}

// Test propagation of periodic orbit from Broucke and Elipe (2005).
// Reference: Broucke, R.A., and Elipe, A. (2005), THE DYNAMICS OF ORBITS IN A POTENTIAL FIELD OF A SOLID CIRCULAR RING,
//              Regular and Chaotic Dynamics, vol 10, num 2, pp 129 - 143
BOOST_AUTO_TEST_CASE( testRingPeriodicOrbit )
{

    // Define ring gravity parameters
    const double gravitationalParameter = 1.0;
    const double ringRadius = 1.0;

    // Define orbit, Broucke and Elipe (2005), section 7
    const Eigen::Vector6d initialState = ( Eigen::Vector6d( ) << 1.8, 0.0, 0.0, 0.0, 0.85873199, 0.0 ).finished();
    const double initialTime = 0.0;
    const double finalTime = 6.58513 * 2.0;

    BodyListSettings bodySettings = BodyListSettings( );

    std::string bodyName = "Ring";
    bodySettings.addSettings( bodyName );
    bodySettings.at( bodyName )->gravityFieldSettings = ringGravitySettings(
                gravitationalParameter, ringRadius, "RingBodyFixed" );
    bodySettings.at( bodyName )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ) );
    bodySettings.at( bodyName )->rotationModelSettings = constantRotationModelSettings(
            "ECLIPJ2000", "RingBodyFixed", Eigen::Matrix3d::Identity() );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 0.0 );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate = { "Vehicle" };
    std::vector< std::string > centralBodies = { "Ring" };

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

    accelerationsOfVehicle[ "Ring" ].push_back( ringAcceleration( ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Define integrator settings.
    const double fixedStepSize = 0.0001;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, initialTime, fixedStepSize );

    // Define propagator settings (Cowell)
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, initialState,
              propagationTimeTerminationSettings( finalTime, true ) );

    // Propagate orbit with Cowell method
    SingleArcDynamicsSimulator< double > dynamicsSimulator( bodies, integratorSettings, propagatorSettings,
                                                            true, false, true );

    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    double computedFinalTime = integrationResult.rbegin( )->first;
    Eigen::Vector6d computedFinalState = integrationResult.rbegin( )->second;

    // For cases in which final state is 0, add 1 (otherwise relative error test will fail)
    Eigen::Vector6d biasVector = ( Eigen::Vector6d( ) << 0.0, 1.0, 1.0, 1.0, 0.0, 1.0 ).finished( );

    Eigen::Vector6d computedFinalStateWithBias = computedFinalState + biasVector;
    Eigen::Vector6d expectedFinalStateWithBias = initialState + biasVector;

    // Error is larger for the y and vx elements (the in-plane elements which should theoretically be zero). For the
    // x and vy elements, tolerance of 1e-6 would be fine.

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedFinalStateWithBias, expectedFinalStateWithBias, 1e-4 );
}

// Test propagation of orbit by comparison with ring described by point masses. Similar to previous test, but here
// a 3-dimensional orbit is used. Additionally, it tests the ring acceleration both for the central and 3rd bodu case.
BOOST_AUTO_TEST_CASE( testRingVersusPointMassesGravityModel )
{

    // Define ring gravity parameters
    const double gravitationalParameter = 1.1;
    const double ringRadius = 0.95;

    const double initialTime = 0.0;
    const Eigen::Vector6d initialState = ( Eigen::Vector6d( ) << 1.8, 0.0, 0.1, 0.0, 0.85873199, 0.0 ).finished();
    const double finalTime = 6.58513 * 0.5;

    // Number of point masses in discrete ring model
    int numPointMasses = 1000;

    for ( std::string ringCentralBody : { "Ring", "COM" } )
    {
        Eigen::Vector6d finalStateRing, finalStatePointMasses;

        for ( unsigned int gravityModelId: { 0, 1 } )
        {
            BodyListSettings bodySettings = BodyListSettings( );

            if ( gravityModelId == 0 )
            {
                std::string bodyName = "Ring";
                bodySettings.addSettings( bodyName );
                bodySettings.at( bodyName )->gravityFieldSettings = ringGravitySettings(
                        gravitationalParameter, ringRadius, "RingBodyFixed", false );
                bodySettings.at( bodyName )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ) );
                bodySettings.at( bodyName )->rotationModelSettings = constantRotationModelSettings(
                        "ECLIPJ2000", "RingBodyFixed", Eigen::Matrix3d::Identity( ) );
            }
            else
            {
                for ( int i = 0; i < numPointMasses; ++i )
                {
                    std::string bodyName = "PointMassBody" + std::to_string( i );
                    bodySettings.addSettings( bodyName );
                    bodySettings.at( bodyName )->gravityFieldSettings = centralGravitySettings(
                            gravitationalParameter / numPointMasses );

                    Eigen::Vector6d state = Eigen::Vector6d::Zero( );
                    state( 0 ) = ringRadius * std::cos( i * 2.0 * mathematical_constants::PI / numPointMasses );
                    state( 1 ) = ringRadius * std::sin( i * 2.0 * mathematical_constants::PI / numPointMasses );

                    bodySettings.at( bodyName )->ephemerisSettings = constantEphemerisSettings( state );
                }

            }

            bodySettings.addSettings( "COM" );
            bodySettings.at( "COM" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ) );

            SystemOfBodies bodies = createSystemOfBodies( bodySettings );

            // Create spacecraft object.
            bodies.createEmptyBody( "Vehicle" );
            bodies.at( "Vehicle" )->setConstantBodyMass( 0.0 );

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate = { "Vehicle" };
            std::vector< std::string > centralBodies;

            // Define propagation settings.
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

            if ( gravityModelId == 0 )
            {
                accelerationsOfVehicle[ "Ring" ].push_back( ringAcceleration( ) );
                centralBodies.push_back( ringCentralBody );
            }
            else
            {
                for ( int i = 0; i < numPointMasses; ++i )
                {
                    std::string bodyName = "PointMassBody" + std::to_string( i );
                    accelerationsOfVehicle[ bodyName ].push_back( pointMassGravityAcceleration( ) );
                }
                centralBodies.push_back( "COM" );
            }
            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

            // Define integrator settings.
            const double fixedStepSize = 0.0001;
            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >( rungeKutta4, initialTime, fixedStepSize );

            // Define propagator settings (Cowell)
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            centralBodies, accelerationModelMap, bodiesToPropagate, initialState,
                            propagationTimeTerminationSettings( finalTime, true ) );

            // Propagate orbit with Cowell method
            SingleArcDynamicsSimulator< double > dynamicsSimulator( bodies, integratorSettings, propagatorSettings,
                                                                    true, false, true );

            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            double computedFinalTime = integrationResult.rbegin( )->first;
            Eigen::Vector6d computedFinalState = integrationResult.rbegin( )->second;

            if ( gravityModelId == 0 )
            {
                finalStateRing = computedFinalState;
            }
            else
            {
                finalStatePointMasses = computedFinalState;
            }

        }

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( finalStateRing, finalStatePointMasses, 1e-14 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace tudat
} // namespace unit_tests