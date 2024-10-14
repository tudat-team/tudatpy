/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Musegaas, P. (2012). Optimization of Space Trajectories Including Multiple Gravity Assists
 *          and Deep Space Maneuvers. MSc Thesis, Delft University of Technology, Delft,
 *          The Netherlands.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <vector>


#include <memory>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <tudat/astro/basic_astro/physicalConstants.h>
#include <tudat/basics/testMacros.h>
#include <tudat/math/basic/mathematicalConstants.h>
#include "tudat/astro/basic_astro/unitConversions.h"
#include <tudat/astro/basic_astro/orbitalElementConversions.h>


#include <tudat/io/basicInputOutput.h>

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/mission_segments/createTransferTrajectory.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include "tudat/math/root_finders/createRootFinder.h"

#include "tudat/astro/low_thrust/shape_based/sphericalShapingLeg.h"
#include "tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h"
#include "tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h"
#include "tudat/astro/low_thrust/shape_based/createBaseFunctionHodographicShaping.h"

namespace tudat
{
namespace unit_tests
{

// Additional namespaces to be used.
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace mission_segments;
using namespace tudat::shape_based_methods;


//! Test implementation of trajectory class
BOOST_AUTO_TEST_SUITE( test_trajectory )

// Create the components for the hodographic shaping radial velocity function used in multiple tests
std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > createHodographicShapingDefaultRadialFunction (
        double scaleFactor, double frequency)
{
    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdRadialVelocityBaseFunctionSettings ) );

    return radialVelocityFunctionComponents;
}

// Create the components for the hodographic shaping axial velocity function used in multiple tests
std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > createHodographicShapingDefaultAxialFunction (
        double scaleFactor, double frequency, int numberOfRevolutions )
{
    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    return axialVelocityFunctionComponents;
}

// Create the components for the hodographic shaping normal velocity function used in multiple tests
std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > createHodographicShapingDefaultNormalFunction (
        double scaleFactor, double frequency)
{
    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );

    return normalVelocityFunctionComponents;
}

BOOST_AUTO_TEST_CASE( testMgaMixedHighLowThrustLegs )
{
    // Set transfer properties
    std::vector< std::string > bodyOrder = { "Earth", "Mars", "Earth", "Venus"};
    int numberOfRevolutions = 2;
    double JD = physical_constants::JULIAN_DAY;
    double departureDate = 6708 * JD;
    std::vector< double > timesOfFlight = { 260.0*JD, 750.0*JD, 300.0*JD };
    double frequency = 2.0 * mathematical_constants::PI / timesOfFlight.at(0);
    double scaleFactor = 1.0 / timesOfFlight.at(0);

    // Create environment
    tudat::simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents =
            createHodographicShapingDefaultRadialFunction ( scaleFactor, frequency );
    // Create components of the axial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents =
            createHodographicShapingDefaultAxialFunction ( scaleFactor, frequency, numberOfRevolutions );
    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents =
            createHodographicShapingDefaultNormalFunction ( scaleFactor, frequency );

    // Create leg and nodes settings
    std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
    std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

    transferLegSettings.resize(bodyOrder.size( ) - 1);
    transferLegSettings[0] = unpoweredLeg( );
    transferLegSettings[1] = hodographicShapingLeg(
                    radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents);
    transferLegSettings[2] = unpoweredLeg( );

    transferNodeSettings.resize(bodyOrder.size( ));
    transferNodeSettings[0] = escapeAndDepartureNode(std::numeric_limits< double >::infinity( ), 0.0);
    transferNodeSettings[1] = swingbyNode();
    transferNodeSettings[2] = swingbyNode();
    transferNodeSettings[3] = captureAndInsertionNode(std::numeric_limits< double >::infinity( ), 0.0);

    // Print parameter definition
    std::cout << "Transfer with unpowered, hodographic shaping, unpowered legs:" << std::endl;
    printTransferParameterDefinition(transferLegSettings, transferNodeSettings);

    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
            bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun");

    // Create list of node times
    std::vector< double > nodeTimes;
    nodeTimes.push_back(departureDate );
    nodeTimes.push_back(nodeTimes.at(0) + timesOfFlight.at(0));
    nodeTimes.push_back(nodeTimes.at(1) + timesOfFlight.at(1));
    nodeTimes.push_back(nodeTimes.at(2) + timesOfFlight.at(2));

    std::vector< Eigen::VectorXd > transferLegFreeParameters(bodyOrder.size( ) - 1);
    transferLegFreeParameters.at(0) = Eigen::VectorXd(0);
    transferLegFreeParameters.at(1) = (Eigen::Vector1d( ) << numberOfRevolutions).finished( );
    transferLegFreeParameters.at(2) = Eigen::VectorXd(0);

    std::vector< Eigen::VectorXd > transferNodeFreeParameters(bodyOrder.size( ));
    // Initial and final excess velocity is 0.0, meaning the spacecraft starts/ends with the velocity of the planet
    double swingbyPeriapsis1 = 6.5e5 * 1e3;
    double swingbyNodeDeltaV1 = 10.0;
    double swingbyRotationAngle1 = 0.0;
    double swingbyPeriapsis2 = 42.0e6;
    double swingbyNodeDeltaV2 = 0.0;
    double swingbyRotationAngle2 = 0.0;

    transferNodeFreeParameters.at(0) = Eigen::VectorXd::Zero( 0 );
    transferNodeFreeParameters.at(1) = ( Eigen::Vector3d( ) << swingbyPeriapsis1, swingbyRotationAngle1, swingbyNodeDeltaV1 ).finished( );
    transferNodeFreeParameters.at(2) = ( Eigen::Vector3d( ) << swingbyPeriapsis2, swingbyRotationAngle2, swingbyNodeDeltaV2 ).finished( );
    transferNodeFreeParameters.at(3) = Eigen::VectorXd::Zero( 0 );

    transferTrajectory->evaluateTrajectory(nodeTimes, transferLegFreeParameters, transferNodeFreeParameters);

    // Check continuity of velocity between legs and nodes
    for( int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
    {
        std::shared_ptr< TransferLeg > leg = transferTrajectory->getLegs().at(i);
        std::shared_ptr< TransferNode > previousNode = transferTrajectory->getNodes().at(i);
        std::shared_ptr< TransferNode > followingNode = transferTrajectory->getNodes().at(i + 1);

        for ( int j = 0; j < 3; j++)
        {
            BOOST_CHECK_SMALL(std::fabs(leg->getDepartureVelocity()[j] - previousNode->getOutgoingVelocity()[j] ), 1e-12);
            BOOST_CHECK_SMALL(std::fabs(leg->getArrivalVelocity()[j] - followingNode->getIncomingVelocity()[j] ), 1e-12);
        }
    }

    // Check whether Delta V at swingby is consistent with node parameters
    BOOST_CHECK_SMALL(std::fabs( transferTrajectory->getNodeDeltaV(1) - swingbyNodeDeltaV1 ), 1e-12);
    BOOST_CHECK_SMALL(std::fabs( transferTrajectory->getNodeDeltaV(2) - swingbyNodeDeltaV2 ), 1e-12);

    // Manually compute outgoing velocity at swingby node and check whether consistent with values in node
    Eigen::Vector3d swingbyNodeVelocity1 = bodies.at( "Mars" )->getEphemeris( )->getCartesianState(
            nodeTimes.at(1) ).segment< 3 >( 3 );
    Eigen::Vector3d swingbyNodeOutgoingVelocity1 = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
            bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( ),
            swingbyNodeVelocity1, transferTrajectory->getNodes().at(1)->getIncomingVelocity(),
            swingbyRotationAngle1, swingbyPeriapsis1, swingbyNodeDeltaV1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( swingbyNodeOutgoingVelocity1, transferTrajectory->getNodes().at(1)->getOutgoingVelocity(), 1.0E-10 );

    Eigen::Vector3d swingbyNodeVelocity2 = bodies.at( "Earth" )->getEphemeris( )->getCartesianState(
            nodeTimes.at(2) ).segment< 3 >( 3 );
    Eigen::Vector3d swingbyNodeIncomingVelocity2 = mission_segments::calculatePoweredGravityAssistIncomingVelocity(
            bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ),
            swingbyNodeVelocity2, transferTrajectory->getNodes().at(2)->getOutgoingVelocity(),
            swingbyRotationAngle2, swingbyPeriapsis2, swingbyNodeDeltaV2 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( swingbyNodeIncomingVelocity2, transferTrajectory->getNodes().at(2)->getIncomingVelocity(), 1.0E-10 );


    // Check if total value of delta V is consistent with manually computed legs
    double totalDeltaV = swingbyNodeDeltaV1 + swingbyNodeDeltaV2;

    for (int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
    {
        if ( i == 1 )
        {
            Eigen::Vector3d departureVelocity = transferTrajectory->getNodes().at(i)->getOutgoingVelocity();
            Eigen::Vector3d arrivalVelocity = transferTrajectory->getNodes().at(i+1)->getIncomingVelocity();
            const std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return departureVelocity; };
            const std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return arrivalVelocity; };

            shape_based_methods::HodographicShapingLeg hodographicShapingLeg = shape_based_methods::HodographicShapingLeg(
                    bodies.at( bodyOrder.at(i) )->getEphemeris( ),
                    bodies.at( bodyOrder.at(i+1) )->getEphemeris( ),
                    bodies.at( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( ),
                    departureVelocityFunction, arrivalVelocityFunction,
                    radialVelocityFunctionComponents, normalVelocityFunctionComponents,
                    axialVelocityFunctionComponents);
            hodographicShapingLeg.updateLegParameters(
                    ( Eigen::Vector3d( )<< nodeTimes.at(i), nodeTimes.at(i+1), numberOfRevolutions ).finished( ) );
            totalDeltaV += hodographicShapingLeg.getLegDeltaV();
        }
        else if ( i == 0 || i == 2)
        {
            mission_segments::UnpoweredUnperturbedTransferLeg unpoweredLeg = mission_segments::UnpoweredUnperturbedTransferLeg(
                    bodies.at( bodyOrder.at(i) )->getEphemeris( ),
                    bodies.at( bodyOrder.at(i+1) )->getEphemeris( ),
                    bodies.at( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( ) );
            unpoweredLeg.updateLegParameters(( Eigen::Vector2d( )<< nodeTimes.at(i), nodeTimes.at(i+1) ).finished( ) );

            if ( i == 0 )
            {
                Eigen::Vector3d departureNodeVelocity = bodies.at(bodyOrder.at(i))->getEphemeris( )->getCartesianState(
                        nodeTimes.at(i)).segment< 3 >(3);
                totalDeltaV += ( unpoweredLeg.getDepartureVelocity( ) - departureNodeVelocity ).norm( );
            }
            else if ( i == 2 )
            {
                Eigen::Vector3d arrivalNodeVelocity = bodies.at(bodyOrder.at(i+1))->getEphemeris( )->getCartesianState(
                        nodeTimes.at(i+1)).segment< 3 >(3);
                totalDeltaV += ( unpoweredLeg.getArrivalVelocity( ) - arrivalNodeVelocity ).norm( );
            }
        }
        else
        {
            throw std::runtime_error( "Unit test with just 3 legs." );
        }

    }
    BOOST_CHECK_SMALL(std::fabs(transferTrajectory->getTotalDeltaV( ) - totalDeltaV ), 1e-12);

    // Retrieve state history and thrust acceleration history.
    transferTrajectory->getStatesAlongTrajectory(50);
    transferTrajectory->getInertialThrustAccelerationsAlongTrajectory(50);

    // Check continuity of retrieved position between legs
    for( int i = 0; i < transferTrajectory->getNumberOfLegs() - 1; i++ )
    {
        std::shared_ptr< TransferLeg > currentLeg = transferTrajectory->getLegs().at(i);
        std::shared_ptr< TransferLeg > nextLeg = transferTrajectory->getLegs().at(i+1);

        Eigen::Vector6d currentLegFinalState = currentLeg->getStateAlongTrajectory( currentLeg->getLegArrivalTime( ) );
        Eigen::Vector6d nextLegInitialState = nextLeg->getStateAlongTrajectory( nextLeg->getLegDepartureTime() );

        // Add 1 to z position because it's value is approximately zero
        Eigen::Vector3d vectorToAdd (0, 0, 1);
        currentLegFinalState.segment(0, 3) += vectorToAdd;
        nextLegInitialState.segment(0, 3) += vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentLegFinalState.segment(0, 3), nextLegInitialState.segment(0, 3), 1e-6);
    }
}

// Test checks delta V of a single hodographic shaping leg, with free parameters for the shaping functions, inserted
// in the MGA framework, comparing it to the results given by Gondelach.
// Also checks for continuity in velocity between leg and nodes
BOOST_AUTO_TEST_CASE( testMgaHodographicShapingSingleLegWithShapingCoefficients )
{

    int numberOfRevolutions = 1;
    double departureDate = 2458849.5; // departureDate is already specified in seconds
    double timeOfFlight = 500.0 * physical_constants::JULIAN_DAY;
    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Set transfer order
    std::vector< std::string > bodyOrder = { "Earth", "Mars" };

    // Create environment
    SystemOfBodies bodies = createSimplifiedSystemOfBodies( );
    // Select correct ephemeris
    bodies.getBody("Earth")->setEphemeris( std::make_shared< ephemerides::ApproximateJplEphemeris >("Earth") );
    bodies.getBody("Mars")->setEphemeris( std::make_shared< ephemerides::ApproximateJplEphemeris >("Mars") );

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, fourthNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, fourthAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, fifthAxialVelocityBaseFunctionSettings ) );


    for( unsigned int creationType = 0; creationType < 2; creationType++ ) {
        // Create leg and nodes settings
        std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
        std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

        if ( creationType == 0 )
        {
            transferLegSettings.resize(bodyOrder.size( ) - 1);
            transferLegSettings[0] = hodographicShapingLeg(
                    radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents);

            transferNodeSettings.resize(bodyOrder.size( ));
            transferNodeSettings[0] = escapeAndDepartureNode(std::numeric_limits< double >::infinity( ), 0.0);
            transferNodeSettings[1] = captureAndInsertionNode(std::numeric_limits< double >::infinity( ), 0.0);

            // Print parameter definition
            std::cout << "Transfer with single hodographic shaping leg with free shaping coefficients:" << std::endl;
            printTransferParameterDefinition(transferLegSettings, transferNodeSettings);
        }
        else if ( creationType == 1 )
        {
            std::vector < shape_based_methods::HodographicBasisFunctionList > radialVelocityFunctionComponentsList;
            radialVelocityFunctionComponentsList.push_back(radialVelocityFunctionComponents);

            std::vector < shape_based_methods::HodographicBasisFunctionList > normalVelocityFunctionComponentsList;
            normalVelocityFunctionComponentsList.push_back(normalVelocityFunctionComponents);

            std::vector < shape_based_methods::HodographicBasisFunctionList > axialVelocityFunctionComponentsList;
            axialVelocityFunctionComponentsList.push_back(axialVelocityFunctionComponents);

            getMgaTransferTrajectorySettingsWithHodographicShapingThrust(
                     transferLegSettings, transferNodeSettings, bodyOrder, radialVelocityFunctionComponentsList,
                     normalVelocityFunctionComponentsList, axialVelocityFunctionComponentsList,
                     std::make_pair(std::numeric_limits< double >::infinity( ), 0.0),
                     std::make_pair(std::numeric_limits< double >::infinity( ), 0.0) );
        }

        std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun");

        // Create list of node times
        std::vector< double > nodeTimes;
        nodeTimes.push_back(departureDate);
        nodeTimes.push_back(nodeTimes.at(0) + timeOfFlight);

        // Create list of leg free parameters
        Eigen::VectorXd freeCoefficientsRadialVelocityFunction = ( Eigen::Vector2d( ) << 500.0, 500.0 ).finished( );
        Eigen::VectorXd freeCoefficientsNormalVelocityFunction = ( Eigen::Vector2d( ) << 500.0, -200.0 ).finished( );
        Eigen::VectorXd freeCoefficientsAxialVelocityFunction = ( Eigen::Vector2d( ) << 500.0, 2000.0 ).finished( );

        std::vector< Eigen::VectorXd > transferLegFreeParameters(bodyOrder.size( ) - 1);
        transferLegFreeParameters.at(0) =
                ( Eigen::Vector7d( ) << numberOfRevolutions, freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction ).finished( );

        // Create list of node free parameters
        std::vector< Eigen::VectorXd > transferNodeFreeParameters(bodyOrder.size( ));
        // Initial and final excess velocity is 0.0, meaning the spacecraft starts/ends with the velocity of the planet
        transferNodeFreeParameters.at(0) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );
        transferNodeFreeParameters.at(1) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );

        transferTrajectory->evaluateTrajectory(nodeTimes, transferLegFreeParameters, transferNodeFreeParameters);

        // Check results consistency w.r.t. results obtained using unit test with original version of Gondelach's code
        double expectedDeltaV = 172313.0;
        BOOST_CHECK_SMALL(std::fabs(transferTrajectory->getTotalDeltaV( ) - expectedDeltaV), 1.0);

        // Check continuity of velocity between legs and nodes
        for( int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
        {
            std::shared_ptr< TransferLeg > leg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferNode > previous_node = transferTrajectory->getNodes().at(i);
            std::shared_ptr< TransferNode > following_node = transferTrajectory->getNodes().at(i+1);

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getDepartureVelocity(), previous_node->getOutgoingVelocity(), 1e-12);
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getArrivalVelocity(), following_node->getIncomingVelocity(), 1e-12);
        }

        // Retrieve state history and thrust acceleration history.
        transferTrajectory->getStatesAlongTrajectory(10);
        transferTrajectory->getInertialThrustAccelerationsAlongTrajectory(10);
    }
}

// Test checks delta V of a single hodographic shaping leg, without free parameters for the shaping functions, inserted
// in the MGA framework, comparing it to the results given by Gondelach.
// Also checks for continuity in velocity between leg and nodes
BOOST_AUTO_TEST_CASE( testMgaHodographicShapingSingleLegWithoutFreeShapingCoefficients )
{

    int numberOfRevolutions = 2;
    double departureDate = 9264.5 * physical_constants::JULIAN_DAY ;
    double timeOfFlight = 1070.0 * physical_constants::JULIAN_DAY;
    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Set transfer order
    std::vector< std::string > bodyOrder = { "Earth", "Mars" };

    // Create environment
    SystemOfBodies bodies = createSimplifiedSystemOfBodies( );
    // Select correct ephemeris
    bodies.getBody("Earth")->setEphemeris( std::make_shared< ephemerides::ApproximateJplEphemeris >("Earth") );
    bodies.getBody("Mars")->setEphemeris( std::make_shared< ephemerides::ApproximateJplEphemeris >("Mars") );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents =
            createHodographicShapingDefaultRadialFunction ( scaleFactor, frequency );
    // Create components of the axial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents =
            createHodographicShapingDefaultAxialFunction ( scaleFactor, frequency, numberOfRevolutions );
    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents =
            createHodographicShapingDefaultNormalFunction ( scaleFactor, frequency );

    for( unsigned int creationType = 0; creationType < 2; creationType++ ) {
        // Create leg and nodes settings
        std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
        std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

        if ( creationType == 0 )
        {
            transferLegSettings.resize(bodyOrder.size( ) - 1);
            transferLegSettings[0] = hodographicShapingLeg(
                    radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents);

            transferNodeSettings.resize(bodyOrder.size( ));
            transferNodeSettings[0] = escapeAndDepartureNode(std::numeric_limits< double >::infinity( ), 0.0);
            transferNodeSettings[1] = captureAndInsertionNode(std::numeric_limits< double >::infinity( ), 0.0);

            // Print parameter definition
            std::cout << "Transfer with single hodographic shaping leg without free shaping coefficients:" << std::endl;
            printTransferParameterDefinition(transferLegSettings, transferNodeSettings);
        }
        else if ( creationType == 1 )
        {
            std::vector < shape_based_methods::HodographicBasisFunctionList > radialVelocityFunctionComponentsList;
            radialVelocityFunctionComponentsList.push_back(radialVelocityFunctionComponents);

            std::vector < shape_based_methods::HodographicBasisFunctionList > normalVelocityFunctionComponentsList;
            normalVelocityFunctionComponentsList.push_back(normalVelocityFunctionComponents);

            std::vector < shape_based_methods::HodographicBasisFunctionList > axialVelocityFunctionComponentsList;
            axialVelocityFunctionComponentsList.push_back(axialVelocityFunctionComponents);

            getMgaTransferTrajectorySettingsWithHodographicShapingThrust(
                     transferLegSettings, transferNodeSettings, bodyOrder, radialVelocityFunctionComponentsList,
                     normalVelocityFunctionComponentsList, axialVelocityFunctionComponentsList,
                     std::make_pair(std::numeric_limits< double >::infinity( ), 0.0),
                     std::make_pair(std::numeric_limits< double >::infinity( ), 0.0) );
        }

        std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun");

        // Create list of node times
        std::vector< double > nodeTimes;
        nodeTimes.push_back(departureDate);
        nodeTimes.push_back(nodeTimes.at(0) + timeOfFlight);

        std::vector< Eigen::VectorXd > transferLegFreeParameters(bodyOrder.size( ) - 1);
        transferLegFreeParameters.at(0) = (Eigen::Vector1d( ) << numberOfRevolutions).finished( );

        std::vector< Eigen::VectorXd > transferNodeFreeParameters(bodyOrder.size( ));
        // Initial and final excess velocity is 0.0, meaning the spacecraft starts/ends with the velocity of the planet
        transferNodeFreeParameters.at(0) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );
        transferNodeFreeParameters.at(1) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );

        transferTrajectory->evaluateTrajectory(nodeTimes, transferLegFreeParameters, transferNodeFreeParameters);

        // Check results consistency w.r.t. Gondelach, D., "A hodographic-shaping method for low-thrust trajectory design",
        // 2012, TU Delft (MSc thesis)
        double expectedDeltaV = 7751.0;
        BOOST_CHECK_SMALL(std::fabs(transferTrajectory->getTotalDeltaV( ) - expectedDeltaV), 1.0);

        // Check continuity of velocity between legs and nodes
        for( int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
        {
            std::shared_ptr< TransferLeg > leg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferNode > previous_node = transferTrajectory->getNodes().at(i);
            std::shared_ptr< TransferNode > following_node = transferTrajectory->getNodes().at(i+1);

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getDepartureVelocity(), previous_node->getOutgoingVelocity(), 1e-12);
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getArrivalVelocity(), following_node->getIncomingVelocity(), 1e-12);
        }

        // Compute peak thrust acceleration
        double peakThrustAcceleration = 0.0;
        std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
        transferTrajectory->getLegs().at(0)->getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

        for( auto it : thrustAccelerationsAlongTrajectory )
        {
            Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
            if ( currentCartesianThrustAcceleration.norm() > peakThrustAcceleration )
            {
                peakThrustAcceleration = currentCartesianThrustAcceleration.norm();
            }
        }

        // Check results consistency w.r.t. Gondelach, D., "A hodographic-shaping method for low-thrust trajectory design",
        // 2012, TU Delft (MSc thesis)
        double expectedPeakAcceleration = 2.64e-4;
        BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-6 );

        // Retrieve state history and thrust acceleration history.
        transferTrajectory->getStatesAlongTrajectory(10);
        transferTrajectory->getInertialThrustAccelerationsAlongTrajectory(10);
    }
}

// Test checks:
// - Continuity of velocity between consecutive nodes and legs
// - Whether deltaV at swingby node is correct (i.e. according to specified as free parameter)
// - Whether outgoing velocity at swingby node is consistent with swinby node incoming velocity and DeltaV
// - Whether the DeltaV computed by manually created hodographic shaping legs is consistent with what is calculated by the transfer trajectory
// Tests assume that the incoming velocity at the swingby node is correct, since it is calculated in the same way as for
// MGA-DSM tests
BOOST_AUTO_TEST_CASE( testMgaMultipleHodographicShapingLegs )
{
    // Set transfer properties
    std::vector< std::string > bodyOrder = { "Earth", "Mars", "Earth" };
    int numberOfRevolutions = 1;
    double JD = physical_constants::JULIAN_DAY;
    double departureDate = 8174.5 * JD;
    std::vector< double > timesOfFlight = { 580.0*JD, 560.0*JD };
    std::vector< double > frequencies = { 2.0 * mathematical_constants::PI / timesOfFlight.at(0), 2.0 * mathematical_constants::PI / timesOfFlight.at(1) };
    std::vector< double > scaleFactors = { 1.0 / timesOfFlight.at(0), 1.0 / timesOfFlight.at(1) };

    // Create environment
    tudat::simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Create components of the radial velocity composite function.
    shape_based_methods::HodographicBasisFunctionList radialVelocityFunctionComponents0 =
            createHodographicShapingDefaultRadialFunction ( scaleFactors.at(0), frequencies.at(0) );
    shape_based_methods::HodographicBasisFunctionList radialVelocityFunctionComponents1 =
            createHodographicShapingDefaultRadialFunction ( scaleFactors.at(1), frequencies.at(1) );
    // Create components of the axial velocity composite function.
    shape_based_methods::HodographicBasisFunctionList axialVelocityFunctionComponents0 =
            createHodographicShapingDefaultAxialFunction ( scaleFactors.at(0), frequencies.at(0), numberOfRevolutions );
    shape_based_methods::HodographicBasisFunctionList axialVelocityFunctionComponents1 =
            createHodographicShapingDefaultAxialFunction ( scaleFactors.at(1), frequencies.at(1), numberOfRevolutions );
    // Create components of the normal velocity composite function.
    shape_based_methods::HodographicBasisFunctionList normalVelocityFunctionComponents0 =
            createHodographicShapingDefaultNormalFunction ( scaleFactors.at(0), frequencies.at(0) );
    shape_based_methods::HodographicBasisFunctionList normalVelocityFunctionComponents1 =
            createHodographicShapingDefaultNormalFunction ( scaleFactors.at(1), frequencies.at(1) );

    for( unsigned int creationType = 0; creationType < 2; creationType++ ) {
        // Create leg and nodes settings
        std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
        std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

        if ( creationType == 0 )
        {
            transferLegSettings.resize(bodyOrder.size( ) - 1);
            transferLegSettings[0] = hodographicShapingLeg(
                    radialVelocityFunctionComponents0, normalVelocityFunctionComponents0, axialVelocityFunctionComponents0);
            transferLegSettings[1] = hodographicShapingLeg(
                    radialVelocityFunctionComponents1, normalVelocityFunctionComponents1, axialVelocityFunctionComponents1);

            transferNodeSettings.resize(bodyOrder.size( ));
            transferNodeSettings[0] = escapeAndDepartureNode(std::numeric_limits< double >::infinity( ), 0.0);
            transferNodeSettings[1] = swingbyNode();
            transferNodeSettings[2] = captureAndInsertionNode(std::numeric_limits< double >::infinity( ), 0.0);

            // Print parameter definition
            std::cout << "Transfer with multiple hodographic shaping legs:" << std::endl;
            printTransferParameterDefinition(transferLegSettings, transferNodeSettings);
        }
        else if ( creationType == 1 )
        {
            std::vector < shape_based_methods::HodographicBasisFunctionList > radialVelocityFunctionComponentsList;
            radialVelocityFunctionComponentsList.push_back(radialVelocityFunctionComponents0);
            radialVelocityFunctionComponentsList.push_back(radialVelocityFunctionComponents1);

            std::vector < shape_based_methods::HodographicBasisFunctionList > normalVelocityFunctionComponentsList;
            normalVelocityFunctionComponentsList.push_back(normalVelocityFunctionComponents0);
            normalVelocityFunctionComponentsList.push_back(normalVelocityFunctionComponents1);

            std::vector < shape_based_methods::HodographicBasisFunctionList > axialVelocityFunctionComponentsList;
            axialVelocityFunctionComponentsList.push_back(axialVelocityFunctionComponents0);
            axialVelocityFunctionComponentsList.push_back(axialVelocityFunctionComponents1);

            getMgaTransferTrajectorySettingsWithHodographicShapingThrust(
                     transferLegSettings, transferNodeSettings, bodyOrder, radialVelocityFunctionComponentsList,
                     normalVelocityFunctionComponentsList, axialVelocityFunctionComponentsList,
                     std::make_pair(std::numeric_limits< double >::infinity( ), 0.0),
                     std::make_pair(std::numeric_limits< double >::infinity( ), 0.0) );
        }

        std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun");

        // Create list of node times
        std::vector< double > nodeTimes;
        nodeTimes.push_back(departureDate);
        nodeTimes.push_back(nodeTimes.at(0) + timesOfFlight.at(0));
        nodeTimes.push_back(nodeTimes.at(1) + timesOfFlight.at(1));

        std::vector< Eigen::VectorXd > transferLegFreeParameters(bodyOrder.size( ) - 1);
        transferLegFreeParameters.at(0) = (Eigen::Vector1d( ) << numberOfRevolutions).finished( );
        transferLegFreeParameters.at(1) = (Eigen::Vector1d( ) << numberOfRevolutions).finished( );

        std::vector< Eigen::VectorXd > transferNodeFreeParameters(bodyOrder.size( ));
        // Initial and final excess velocity is 0.0, meaning the spacecraft starts/ends with the velocity of the planet
        double swingbyNodeDeltaV = 10.0;
        double swingbyRotationAngle = 1.3;
        double swingbyPeriapsis = 65000.0e3;
        transferNodeFreeParameters.at(0) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );
        transferNodeFreeParameters.at(1) = ( Eigen::Vector6d( ) << 10, 0.0, 0.0, swingbyPeriapsis, swingbyRotationAngle, swingbyNodeDeltaV ).finished( );
        transferNodeFreeParameters.at(2) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );

        transferTrajectory->evaluateTrajectory(nodeTimes, transferLegFreeParameters, transferNodeFreeParameters);

        // Check continuity of velocity between legs and nodes
        for( int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
        {
            std::shared_ptr< TransferLeg > leg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferNode > previousNode = transferTrajectory->getNodes().at(i);
            std::shared_ptr< TransferNode > followingNode = transferTrajectory->getNodes().at(i + 1);

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getDepartureVelocity(), previousNode->getOutgoingVelocity(), 1e-12);
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getArrivalVelocity(), followingNode->getIncomingVelocity(), 1e-12);
        }


        // Check whether Delta V at swingby is consistent with node parameters
        BOOST_CHECK_SMALL(std::fabs( transferTrajectory->getNodeDeltaV(1) - swingbyNodeDeltaV ), 1e-12);

        // Manually compute outgoing velocity at swingby node and check whether consistent with values in node
        Eigen::Vector3d swingbyNodeVelocity = bodies.at("Mars" )->getEphemeris( )->getCartesianState(
                nodeTimes.at(0) + timesOfFlight.at(0) ).segment< 3 >( 3 );
        Eigen::Vector3d swingbyNodeOutgoingVelocity = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( ),
                swingbyNodeVelocity, transferTrajectory->getNodes().at(1)->getIncomingVelocity(),
                swingbyRotationAngle, swingbyPeriapsis, swingbyNodeDeltaV );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(swingbyNodeOutgoingVelocity, transferTrajectory->getNodes().at(1)->getOutgoingVelocity(), 1e-12);


        // Check if total value of delta V is correct
        double totalDeltaV = swingbyNodeDeltaV;

        for (int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
        {
            Eigen::Vector3d departureVelocity;
            Eigen::Vector3d arrivalVelocity;
            if (i == 0)
            {
                departureVelocity = bodies.at(bodyOrder.at(i) )->getEphemeris( )->getCartesianState(
                        nodeTimes.at(i) ).segment< 3 >( 3 );
                arrivalVelocity = transferTrajectory->getNodes().at(i+1)->getIncomingVelocity();
            }
            else if (i == 1)
            {
                departureVelocity = transferTrajectory->getNodes().at(i)->getOutgoingVelocity();
                arrivalVelocity = bodies.at(bodyOrder.at(i+1) )->getEphemeris( )->getCartesianState(
                        nodeTimes.at(i+1) ).segment< 3 >( 3 );
            }
            else
            {
                throw std::runtime_error( "Unit test with just 2 legs." );
            }

            const std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return departureVelocity; };
            const std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return arrivalVelocity; };

            shape_based_methods::HodographicBasisFunctionList legRadialVelocityFunctionComponents;
            shape_based_methods::HodographicBasisFunctionList legNormalVelocityFunctionComponents;
            shape_based_methods::HodographicBasisFunctionList legAxialVelocityFunctionComponents;

            if ( i == 0 )
            {
                legRadialVelocityFunctionComponents = radialVelocityFunctionComponents0;
                legNormalVelocityFunctionComponents = normalVelocityFunctionComponents0;
                legAxialVelocityFunctionComponents = axialVelocityFunctionComponents0;
            }
            else if ( i == 1 )
            {
                legRadialVelocityFunctionComponents = radialVelocityFunctionComponents1;
                legNormalVelocityFunctionComponents = normalVelocityFunctionComponents1;
                legAxialVelocityFunctionComponents = axialVelocityFunctionComponents1;
            }

            shape_based_methods::HodographicShapingLeg hodographicShapingLeg = shape_based_methods::HodographicShapingLeg(
                    bodies.at( bodyOrder.at(i) )->getEphemeris( ),
                    bodies.at( bodyOrder.at(i+1) )->getEphemeris( ),
                    bodies.at( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( ),
                    departureVelocityFunction, arrivalVelocityFunction,
                    legRadialVelocityFunctionComponents, legNormalVelocityFunctionComponents,
                    legAxialVelocityFunctionComponents);


            hodographicShapingLeg.updateLegParameters(
                    ( Eigen::Vector3d( )<< nodeTimes.at(i), nodeTimes.at(i+1), numberOfRevolutions ).finished( ) );

            totalDeltaV += hodographicShapingLeg.getLegDeltaV();
        }
        BOOST_CHECK_SMALL(std::fabs(transferTrajectory->getTotalDeltaV( ) - totalDeltaV ), 1e-12);

        // Retrieve state history and thrust acceleration history.
        transferTrajectory->getStatesAlongTrajectory(10);
        transferTrajectory->getInertialThrustAccelerationsAlongTrajectory(10);

        // Check continuity of retrieved position between legs
        for( int i = 0; i < transferTrajectory->getNumberOfLegs() - 1; i++ )
        {
            std::shared_ptr< TransferLeg > currentLeg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferLeg > nextLeg = transferTrajectory->getLegs().at(i+1);

            Eigen::Vector6d currentLegFinalState = currentLeg->getStateAlongTrajectory( currentLeg->getLegArrivalTime( ) );
            Eigen::Vector6d nextLegInitialState = nextLeg->getStateAlongTrajectory( nextLeg->getLegDepartureTime() );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentLegFinalState.segment(0, 3), nextLegInitialState.segment(0, 3), 1e-12);

        }
    }
}

// Test checks delta V of a single spherical shaping leg, inserted in the MGA framework, comparing it to the results
// given by Roegiers.
// Also checks for continuity in velocity between leg and nodes
BOOST_AUTO_TEST_CASE( testMgaSphericalShapingSingleLeg )
{
    int numberOfRevolutions = 1;
    double departureDate = 8174.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 580.0 * physical_constants::JULIAN_DAY;

    // Create environment
    tudat::simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Set transfer order
    std::vector< std::string > bodyOrder = { "Earth", "Mars" };

    // Define root finder settings
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            tudat::root_finders::bisectionRootFinderSettings( 1.0E-6, TUDAT_NAN, TUDAT_NAN, 30 );
    double lowerBoundFreeCoefficient = 1.0e-6;
    double upperBoundFreeCoefficient = 1.0e-1;

    for( unsigned int creationType = 0; creationType < 2; creationType++ ) {
        // Create leg and nodes settings
        std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
        std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

        if ( creationType == 0 )
        {
            transferLegSettings.resize(bodyOrder.size( ) - 1);
            transferLegSettings[0] = sphericalShapingLeg(rootFinderSettings,
                                                         lowerBoundFreeCoefficient, upperBoundFreeCoefficient);

            transferNodeSettings.resize(bodyOrder.size( ));
            transferNodeSettings[0] = escapeAndDepartureNode(std::numeric_limits< double >::infinity( ), 0.0);
            transferNodeSettings[1] = captureAndInsertionNode(std::numeric_limits< double >::infinity( ), 0.0);

            // Print parameter definition
            std::cout << "Transfer with single spherical shaping leg:" << std::endl;
            printTransferParameterDefinition(transferLegSettings, transferNodeSettings);
        }
        else if ( creationType == 1 )
        {
            getMgaTransferTrajectorySettingsWithSphericalShapingThrust(
                    transferLegSettings, transferNodeSettings, bodyOrder, rootFinderSettings,
                    std::make_pair(std::numeric_limits< double >::infinity( ), 0.0),
                    std::make_pair(std::numeric_limits< double >::infinity( ), 0.0),
                    lowerBoundFreeCoefficient, upperBoundFreeCoefficient);
        }

        std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun");

        // Create list of node times
        std::vector< double > nodeTimes;
        nodeTimes.push_back(departureDate);
        nodeTimes.push_back(nodeTimes.at(0) + timeOfFlight);

        std::vector< Eigen::VectorXd > transferLegFreeParameters(bodyOrder.size( ) - 1);
        transferLegFreeParameters.at(0) = (Eigen::Vector1d( ) << numberOfRevolutions).finished( );

        std::vector< Eigen::VectorXd > transferNodeFreeParameters(bodyOrder.size( ));
        // Initial and final excess velocity is 0.0, meaning the spacecraft starts/ends with the velocity of the planet
        transferNodeFreeParameters.at(0) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );
        transferNodeFreeParameters.at(1) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );

        transferTrajectory->evaluateTrajectory(nodeTimes, transferLegFreeParameters, transferNodeFreeParameters);

        // Check results consistency w.r.t. Roegiers, T., Application of the Spherical Shaping Method to a Low-Thrust
        // Multiple Asteroid Rendezvous Mission, TU Delft (MSc thesis), 2014
        double expectedDeltaV = 5700.0;
        BOOST_CHECK_SMALL(std::fabs(transferTrajectory->getTotalDeltaV( ) - expectedDeltaV), 5.0);

        // Check continuity of velocity between legs and nodes
        for( int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
        {
            std::shared_ptr< TransferLeg > leg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferNode > previous_node = transferTrajectory->getNodes().at(i);
            std::shared_ptr< TransferNode > following_node = transferTrajectory->getNodes().at(i+1);

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getDepartureVelocity(), previous_node->getOutgoingVelocity(), 1e-12);
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getArrivalVelocity(), following_node->getIncomingVelocity(), 1e-12);
        }

        // Retrieve state history and thrust acceleration history.
        transferTrajectory->getStatesAlongTrajectory(10);
        transferTrajectory->getInertialThrustAccelerationsAlongTrajectory(10);
        transferTrajectory->getRswThrustAccelerationsAlongTrajectory(10);
        std::map< double, Eigen::Vector3d > tnwAcceleration = transferTrajectory->getTnwThrustAccelerationsAlongTrajectory(50);

        for( std::map< double, Eigen::Vector3d >::iterator it = tnwAcceleration.begin(); it != tnwAcceleration.end(); ++it )
        {
            // Check if thrust in N direction is zero (it should be, see Roegiers, section 6.5)
            BOOST_CHECK_SMALL(std::fabs(it->second.transpose()[1]), 1e-16);
        }
    }
}

// Test checks:
// - Continuity of velocity between consecutive nodes and legs
// - Whether deltaV at swingby node is correct (i.e. according to specified as free parameter)
// - Whether outgoing velocity at swingby node is consistent with swinby node incoming velocity and DeltaV
// - Whether the DeltaV computed by manually created spherical shaping legs is consistent with what is calculated by the transfer trajectory
// Tests assume that the incoming velocity at the swingby node is correct, since it is calculated in the same way as for
// MGA-DSM tests
BOOST_AUTO_TEST_CASE( testMgaMultipleSphericalShapingLegs )
{
    // Set transfer properties
    std::vector< std::string > bodyOrder = { "Earth", "Mars", "Earth" };
    int numberOfRevolutions = 1;
    double JD = physical_constants::JULIAN_DAY;
    double departureDate = 8174.5 * JD;
    std::vector< double > timesOfFlight = { 580.0*JD, 560.0*JD };

    // Create environment
    tudat::simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Define root finder settings
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            tudat::root_finders::bisectionRootFinderSettings( 1.0E-6, TUDAT_NAN, TUDAT_NAN, 30 );
    double lowerBoundFreeCoefficient = 1.0e-6;
    double upperBoundFreeCoefficient = 1.0e-1;

    for( unsigned int creationType = 0; creationType < 2; creationType++ ) {
        // Create leg and nodes settings
        std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
        std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

        if ( creationType == 0 )
        {
            transferLegSettings.resize(bodyOrder.size( ) - 1);
            transferLegSettings[0] = sphericalShapingLeg(rootFinderSettings,
                                                         lowerBoundFreeCoefficient, upperBoundFreeCoefficient);
            transferLegSettings[1] = sphericalShapingLeg(rootFinderSettings,
                                                         lowerBoundFreeCoefficient, upperBoundFreeCoefficient);

            transferNodeSettings.resize(bodyOrder.size( ));
            transferNodeSettings[0] = escapeAndDepartureNode(std::numeric_limits< double >::infinity( ), 0.0);
            transferNodeSettings[1] = swingbyNode();
            transferNodeSettings[2] = captureAndInsertionNode(std::numeric_limits< double >::infinity( ), 0.0);

            // Print parameter definition
            std::cout << "Transfer with multiple spherical shaping legs:" << std::endl;
            printTransferParameterDefinition(transferLegSettings, transferNodeSettings);
        }
        else if ( creationType == 1 )
        {
            getMgaTransferTrajectorySettingsWithSphericalShapingThrust(
                    transferLegSettings, transferNodeSettings, bodyOrder, rootFinderSettings,
                    std::make_pair(std::numeric_limits< double >::infinity( ), 0.0),
                    std::make_pair(std::numeric_limits< double >::infinity( ), 0.0),
                    lowerBoundFreeCoefficient, upperBoundFreeCoefficient);
        }

        std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun");

        // Create list of node times
        std::vector< double > nodeTimes;
        nodeTimes.push_back(departureDate);
        nodeTimes.push_back(nodeTimes.at(0) + timesOfFlight.at(0));
        nodeTimes.push_back(nodeTimes.at(1) + timesOfFlight.at(1));

        std::vector< Eigen::VectorXd > transferLegFreeParameters(bodyOrder.size( ) - 1);
        transferLegFreeParameters.at(0) = (Eigen::Vector1d( ) << numberOfRevolutions).finished( );
        transferLegFreeParameters.at(1) = (Eigen::Vector1d( ) << numberOfRevolutions).finished( );

        std::vector< Eigen::VectorXd > transferNodeFreeParameters(bodyOrder.size( ));
        // Initial and final excess velocity is 0.0, meaning the spacecraft starts/ends with the velocity of the planet
        double swingbyNodeDeltaV = 10.0;
        double swingbyRotationAngle = 1.3;
        double swingbyPeriapsis = 65000.0e3;
        transferNodeFreeParameters.at(0) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );
        transferNodeFreeParameters.at(1) = ( Eigen::Vector6d( ) << 10, 0.0, 0.0, swingbyPeriapsis, swingbyRotationAngle, swingbyNodeDeltaV ).finished( );
        transferNodeFreeParameters.at(2) = ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( );

        transferTrajectory->evaluateTrajectory(nodeTimes, transferLegFreeParameters, transferNodeFreeParameters);

        // Check continuity of velocity between legs and nodes
        for( int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
        {
            std::shared_ptr< TransferLeg > leg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferNode > previousNode = transferTrajectory->getNodes().at(i);
            std::shared_ptr< TransferNode > followingNode = transferTrajectory->getNodes().at(i + 1);

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getDepartureVelocity(), previousNode->getOutgoingVelocity(), 1e-12);
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(leg->getArrivalVelocity(), followingNode->getIncomingVelocity(), 1e-12);
        }


        // Check whether Delta V at swingby is consistent with node parameters
        BOOST_CHECK_SMALL(std::fabs( transferTrajectory->getNodeDeltaV(1) - swingbyNodeDeltaV ), 1e-12);

        // Manually compute outgoing velocity at swingby node and check whether consistent with values in node
        Eigen::Vector3d swingbyNodeVelocity = bodies.at("Mars" )->getEphemeris( )->getCartesianState(
                nodeTimes.at(0) + timesOfFlight.at(0) ).segment< 3 >( 3 );
        Eigen::Vector3d swingbyNodeOutgoingVelocity = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( ),
                swingbyNodeVelocity, transferTrajectory->getNodes().at(1)->getIncomingVelocity(),
                swingbyRotationAngle, swingbyPeriapsis, swingbyNodeDeltaV );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(swingbyNodeOutgoingVelocity, transferTrajectory->getNodes().at(1)->getOutgoingVelocity(), 1e-12);

        // Check if total value of delta V is correct
        double totalDeltaV = swingbyNodeDeltaV;

        for (int i = 0; i < transferTrajectory->getNumberOfLegs(); i++ )
        {
            Eigen::Vector3d departureVelocity;
            Eigen::Vector3d arrivalVelocity;
            if (i == 0)
            {
                departureVelocity = bodies.at(bodyOrder.at(i) )->getEphemeris( )->getCartesianState(
                        nodeTimes.at(i) ).segment< 3 >( 3 );
                arrivalVelocity = transferTrajectory->getNodes().at(i+1)->getIncomingVelocity();
            }
            else if (i == 1)
            {
                departureVelocity = transferTrajectory->getNodes().at(i)->getOutgoingVelocity();
                arrivalVelocity = bodies.at(bodyOrder.at(i+1) )->getEphemeris( )->getCartesianState(
                        nodeTimes.at(i+1) ).segment< 3 >( 3 );
            }
            else
            {
                throw std::runtime_error( "Unit test with just 2 legs." );
            }

            const std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return departureVelocity; };
            const std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return arrivalVelocity; };

            shape_based_methods::SphericalShapingLeg sphericalShapingLeg = shape_based_methods::SphericalShapingLeg (
                    bodies.at( bodyOrder.at(i) )->getEphemeris( ),
                    bodies.at( bodyOrder.at(i+1) )->getEphemeris( ),
                    bodies.at( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( ),
                    departureVelocityFunction, arrivalVelocityFunction,
                    rootFinderSettings, lowerBoundFreeCoefficient, upperBoundFreeCoefficient);

            sphericalShapingLeg.updateLegParameters(
                    ( Eigen::Vector3d( )<< nodeTimes.at(i), nodeTimes.at(i+1), numberOfRevolutions ).finished( ) );

            totalDeltaV += sphericalShapingLeg.getLegDeltaV();
        }
        BOOST_CHECK_SMALL(std::fabs(transferTrajectory->getTotalDeltaV( ) - totalDeltaV ), 1e-12);

        // Retrieve state history and thrust acceleration history.
        transferTrajectory->getStatesAlongTrajectory(10);
        transferTrajectory->getInertialThrustAccelerationsAlongTrajectory(10);
        transferTrajectory->getRswThrustAccelerationsAlongTrajectory(10);
        std::map< double, Eigen::Vector3d > tnwAcceleration = transferTrajectory->getTnwThrustAccelerationsAlongTrajectory(50);

        for( std::map< double, Eigen::Vector3d >::iterator it = tnwAcceleration.begin(); it != tnwAcceleration.end(); ++it )
        {
            // Check if thrust in N direction is zero (it should be, see Roegiers, section 6.5)
            BOOST_CHECK_SMALL(std::fabs(it->second.transpose()[1]), 1e-16);
        }

        // Check continuity of retrieved position between legs
        for( int i = 0; i < transferTrajectory->getNumberOfLegs() - 1; i++ )
        {
            std::shared_ptr< TransferLeg > currentLeg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferLeg > nextLeg = transferTrajectory->getLegs().at(i+1);

            Eigen::Vector6d currentLegFinalState = currentLeg->getStateAlongTrajectory( currentLeg->getLegArrivalTime( ) );
            Eigen::Vector6d nextLegInitialState = nextLeg->getStateAlongTrajectory( nextLeg->getLegDepartureTime() );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentLegFinalState.segment(0, 3), nextLegInitialState.segment(0, 3), 1e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE( testMGATrajectory_New )
{
    // Expected test result based on the ideal Cassini 1 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 4930.72686847243;

    // Create environment
    tudat::simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Set transfer order
    std::vector< std::string > bodyOrder = {
        "Earth", "Venus", "Venus", "Earth", "Jupiter", "Saturn" };
    int numberOfNodes = bodyOrder.size( );

    std::shared_ptr< TransferTrajectory > transferTrajectory;
    double nominalDeltaV = TUDAT_NAN;
    double nominalCaptureDeltaV = TUDAT_NAN;

    for( unsigned int creationType = 0; creationType < 3; creationType++ )
    {
        // Create leg settings (all unpowered)
        std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
        std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;

        // Define minimum periapsis altitudes for flybys;
        std::map< std::string, double > minimumPeriapses;
        minimumPeriapses[ "Venus" ] = 6351800.0;
        minimumPeriapses[ "Earth" ] = 6678000.0;
        minimumPeriapses[ "Jupiter" ] =  600000000.0;
        minimumPeriapses[ "Saturn" ] = 65000000.0;

        if( creationType == 0 )
        {
            transferLegSettings.resize( numberOfNodes - 1 );
            transferLegSettings[ 0 ] = unpoweredLeg( );
            transferLegSettings[ 1 ] = unpoweredLeg( );
            transferLegSettings[ 2 ] = unpoweredLeg( );
            transferLegSettings[ 3 ] = unpoweredLeg( );
            transferLegSettings[ 4 ] = unpoweredLeg( );

            transferNodeSettings.resize( numberOfNodes );
            transferNodeSettings[ 0 ] = escapeAndDepartureNode( std::numeric_limits< double >::infinity( ), 0.0 );
            transferNodeSettings[ 1 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 1 ) ) );
            transferNodeSettings[ 2 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 2 ) ) );
            transferNodeSettings[ 3 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 3 ) ) );
            transferNodeSettings[ 4 ] = swingbyNode( minimumPeriapses.at( bodyOrder.at( 4 ) ) );
            transferNodeSettings[ 5 ] = captureAndInsertionNode( 1.0895e8 / 0.02, 0.98 );
        }
        else if( creationType == 1 )
        {
            getMgaTransferTrajectorySettingsWithoutDsm(
                        transferLegSettings, transferNodeSettings, bodyOrder,
                        std::make_pair( std::numeric_limits< double >::infinity( ), 0.0 ),
                        std::make_pair( 1.0895e8 / 0.02, 0.98 ),
                        minimumPeriapses );
        }
        else if( creationType == 2 )
        {
            getMgaTransferTrajectorySettingsWithoutDsm(
                        transferLegSettings, transferNodeSettings, bodyOrder,
                        std::make_pair( std::numeric_limits< double >::infinity( ), 0.0 ),
                        std::make_pair( TUDAT_NAN, TUDAT_NAN ),
                        minimumPeriapses );

        }
        transferTrajectory = createTransferTrajectory(
                    bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun" );




        // Create list of node times
        double JD = physical_constants::JULIAN_DAY;
        std::vector< double > nodeTimes;
        // Subtract 0.5 to go from MJD2000 to J2000
        nodeTimes.push_back( ( -789.8117 - 0.5 ) * JD );
        nodeTimes.push_back( nodeTimes.at( 0 ) + 158.302027105278 * JD );
        nodeTimes.push_back( nodeTimes.at( 1 ) + 449.385873819743 * JD );
        nodeTimes.push_back( nodeTimes.at( 2 ) + 54.7489684339665 * JD );
        nodeTimes.push_back( nodeTimes.at( 3 ) + 1024.36205846918 * JD );
        nodeTimes.push_back( nodeTimes.at( 4 ) + 4552.30796805542 * JD );

        std::vector< Eigen::VectorXd > transferLegFreeParameters;
        for( int i = 0; i < numberOfNodes - 1; i++ )
        {
            transferLegFreeParameters.push_back( Eigen::VectorXd( 0 ) );
        }
        std::vector< Eigen::VectorXd > transferNodeFreeParameters;
        for( int i = 0; i < numberOfNodes; i++ )
        {
            transferNodeFreeParameters.push_back( Eigen::VectorXd( 0 ) );
        }
        if( creationType == 2 )
        {
            transferNodeFreeParameters[ numberOfNodes - 1 ] = ( Eigen::Vector3d( )<<65000000.0, 0.0, 0.0 ).finished( );
        }

        std::cout << "Transfer with multiple unpowered legs:" << std::endl;
        printTransferParameterDefinition( transferLegSettings, transferNodeSettings );

        transferTrajectory->evaluateTrajectory(
                    nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );

        if( creationType < 2 )
        {
            //BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-3 );
            BOOST_CHECK( std::fabs( expectedDeltaV - transferTrajectory->getTotalDeltaV( ) ) < std::fabs( std::min( expectedDeltaV,transferTrajectory->getTotalDeltaV( ) ) * 1.0E-3 ) );

            if( creationType == 0 )
            {
                nominalDeltaV = transferTrajectory->getTotalDeltaV( );
                nominalCaptureDeltaV = transferTrajectory->getNodeDeltaV( 5 );
            }
            else
            {
                //BOOST_CHECK_CLOSE_FRACTION( nominalDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-12 );
                BOOST_CHECK( std::fabs( nominalDeltaV - transferTrajectory->getTotalDeltaV( ) ) < std::fabs( std::min( nominalDeltaV,transferTrajectory->getTotalDeltaV( ) ) * 1.0E-12 ) );

            }
        }
        else if( creationType == 2 )
        {

            //BOOST_CHECK_CLOSE_FRACTION( nominalCaptureDeltaV, nominalDeltaV - transferTrajectory->getTotalDeltaV( ), 1.0E-12 );
            BOOST_CHECK( std::fabs( nominalCaptureDeltaV - ( nominalDeltaV - transferTrajectory->getTotalDeltaV( ) ) ) < std::fabs( std::min( nominalCaptureDeltaV,( nominalDeltaV - transferTrajectory->getTotalDeltaV( ) ) ) * 1.0E-12 ) );
        }

        std::vector< std::map< double, Eigen::Vector6d > > statesAlongTrajectoryPerLeg;
        transferTrajectory->getStatesAlongTrajectoryPerLeg(
                    statesAlongTrajectoryPerLeg, 10 );
        std::vector< std::map< double, Eigen::Vector3d > > thrustAccelerationsAlongTrajectoryPerLeg;
        transferTrajectory->getInertialThrustAccelerationsAlongTrajectoryPerLeg(
                thrustAccelerationsAlongTrajectoryPerLeg, 10);

        double sunGravitationalParameter = bodies.at( "Sun" )->getGravitationalParameter( );
        for( unsigned int i = 0; i < statesAlongTrajectoryPerLeg.size( ); i++ )
        {
            double legStartTime = nodeTimes.at( i );
            double legEndTime = nodeTimes.at( i + 1 );

            std::map< double, Eigen::Vector6d > statesAlongSingleLeg = statesAlongTrajectoryPerLeg.at( i );

            // Check initial and final time on output list
            //BOOST_CHECK_CLOSE_FRACTION( statesAlongSingleLeg.begin( )->first, legStartTime, 1.0E-14 );
            //BOOST_CHECK_CLOSE_FRACTION( statesAlongSingleLeg.rbegin( )->first, legEndTime, 1.0E-14 )
            BOOST_CHECK( std::fabs( statesAlongSingleLeg.begin( )->first - legStartTime ) < std::fabs( std::min( statesAlongSingleLeg.begin( )->first,legStartTime ) * 1.0E-14 ) );
            BOOST_CHECK( std::fabs( statesAlongSingleLeg.rbegin( )->first - legEndTime ) < std::fabs( std::min( statesAlongSingleLeg.rbegin( )->first,legEndTime ) * 1.0E-14 ) );


            // Check if Keplerian state (slow elements) is the same for each output point
            Eigen::Vector6d previousKeplerianState = Eigen::Vector6d::Constant( TUDAT_NAN );
            for( auto it : statesAlongSingleLeg )
            {
                Eigen::Vector6d currentCartesianState = it.second;
                Eigen::Vector6d currentKeplerianState = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                            currentCartesianState, sunGravitationalParameter );
                if( previousKeplerianState == previousKeplerianState )
                {
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                ( previousKeplerianState.segment( 0, 5 ) ),
                                ( currentKeplerianState.segment( 0, 5 ) ),
                                1.0E-10 );

                }
                previousKeplerianState = currentKeplerianState;
            }

            // Check if output meets boundary conditions
            Eigen::Vector3d depatureBodyPosition = bodies.at( bodyOrder.at( i ) )->getStateInBaseFrameFromEphemeris(
                        legStartTime ).segment( 0, 3 );
            Eigen::Vector3d arrivalBodyPosition = bodies.at( bodyOrder.at( i + 1 ) )->getStateInBaseFrameFromEphemeris(
                        legEndTime ).segment( 0, 3 );

            for( int j = 0; j < 3; j++ )
            {
                //TODO: Find out why tolerance needs to be so big for one of the legs
                //BOOST_CHECK_SMALL( std::fabs( statesAlongSingleLeg.begin( )->second( j ) - depatureBodyPosition( j ) ), 20.0E3 );
                //BOOST_CHECK_SMALL( std::fabs( statesAlongSingleLeg.rbegin( )->second( j ) - arrivalBodyPosition( j ) ), 20.0E3 );
                BOOST_CHECK_LT( std::fabs( statesAlongSingleLeg.begin( )->second( j ) - depatureBodyPosition( j ) ), 20.0E3 );
                BOOST_CHECK_LT( std::fabs( statesAlongSingleLeg.rbegin( )->second( j ) - arrivalBodyPosition( j ) ), 30.0E3 );
            }
        }

        // Check continuity of retrieved position between legs
        for( int i = 0; i < transferTrajectory->getNumberOfLegs() - 1; i++ )
        {
            std::shared_ptr< TransferLeg > currentLeg = transferTrajectory->getLegs().at(i);
            std::shared_ptr< TransferLeg > nextLeg = transferTrajectory->getLegs().at(i+1);

            Eigen::Vector6d currentLegFinalState = currentLeg->getStateAlongTrajectory( currentLeg->getLegArrivalTime( ) );
            Eigen::Vector6d nextLegInitialState = nextLeg->getStateAlongTrajectory( nextLeg->getLegDepartureTime() );

            // Add 1 to z position because it's value is approximately zero
            Eigen::Vector3d vectorToAdd (0, 0, 1);
            currentLegFinalState.segment(0, 3) += vectorToAdd;
            nextLegInitialState.segment(0, 3) += vectorToAdd;
            // Not sure why tolerance needs to be so large... probably related to why the tolerance above is also large
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentLegFinalState.segment(0, 3), nextLegInitialState.segment(0, 3), 1e-5);
        }
    }

}


//! Test delta-V computation for MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory1 )
{
    // Expected test result based on the ideal Messenger trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8630.83256199051;

    // Create environment
    simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Set transfer order
    std::vector< std::string > bodyOrder = {
        "Earth", "Earth", "Venus", "Venus", "Mercury" };
    int numberOfNodes = bodyOrder.size( );

    // Create leg settings
    std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
    transferLegSettings.resize( numberOfNodes - 1 );
    transferLegSettings[ 0 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 1 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 2 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 3 ] = dsmVelocityBasedLeg( );

    // Create node settings
    std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;
    transferNodeSettings.resize( numberOfNodes );
    transferNodeSettings[ 0 ] = escapeAndDepartureNode( std::numeric_limits< double >::infinity( ), 0.0 );
    transferNodeSettings[ 1 ] = swingbyNode( );
    transferNodeSettings[ 2 ] = swingbyNode( );
    transferNodeSettings[ 3 ] = swingbyNode( );
    transferNodeSettings[ 4 ] = captureAndInsertionNode( std::numeric_limits< double >::infinity( ), 0.0 );

    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun" );

    // Add the time of flight and start epoch, which are in JD.
    double JD = physical_constants::JULIAN_DAY;
    std::vector< double > nodeTimes;
    nodeTimes.push_back( ( 1171.64503236 - 0.5 ) * JD );
    nodeTimes.push_back( nodeTimes.at( 0 ) + 399.999999715 * JD );
    nodeTimes.push_back( nodeTimes.at( 1 ) + 178.372255301 * JD );
    nodeTimes.push_back( nodeTimes.at( 2 ) + 299.223139512  * JD );
    nodeTimes.push_back( nodeTimes.at( 3 ) + 180.510754824 * JD );

    std::vector< Eigen::VectorXd > transferLegFreeParameters;
    transferLegFreeParameters.resize( numberOfNodes - 1 );
    transferLegFreeParameters[ 0 ] = ( Eigen::VectorXd( 1 ) << 0.234594654679 ).finished( );
    transferLegFreeParameters[ 1 ] = ( Eigen::VectorXd( 1 ) << 0.0964769387134 ).finished( );
    transferLegFreeParameters[ 2 ] = ( Eigen::VectorXd( 1 ) << 0.829948744508 ).finished( );
    transferLegFreeParameters[ 3 ] = ( Eigen::VectorXd( 1 ) << 0.317174785637 ).finished( );

    std::vector< Eigen::VectorXd > transferNodeFreeParameters;
    transferNodeFreeParameters.resize( numberOfNodes );
    transferNodeFreeParameters[ 0 ] = ( Eigen::VectorXd( 3 ) <<1408.99421278, 0.37992647165 * 2 * 3.14159265358979, std::acos(  2 * 0.498004040298 - 1. ) - 3.14159265358979 / 2 ).finished( );
    transferNodeFreeParameters[ 1 ] = ( Eigen::VectorXd( 3 ) <<1.80629232251 * 6.378e6, 1.35077257078, 0.0 ).finished( );
    transferNodeFreeParameters[ 2 ] = ( Eigen::VectorXd( 3 ) <<3.04129845698 * 6.052e6, 1.09554368115, 0.0 ).finished( );
    transferNodeFreeParameters[ 3 ] = ( Eigen::VectorXd( 3 ) <<1.10000000891 * 6.052e6, 1.34317576594, 0.0 ).finished( );


    transferTrajectory->evaluateTrajectory(
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    std::cout << "Transfer with multiple velocity-based DSM legs:" << std::endl;
    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );

    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-3 );

    // Check continuity of retrieved position between legs
    // Tolerance is large because the state is retrieved from the converged state history (so not from the leg's
    // boundary conditions)
    for( int i = 0; i < transferTrajectory->getNumberOfLegs() - 1; i++ )
    {
        std::shared_ptr< TransferLeg > currentLeg = transferTrajectory->getLegs().at(i);
        std::shared_ptr< TransferLeg > nextLeg = transferTrajectory->getLegs().at(i+1);

        Eigen::Vector6d currentLegFinalState = currentLeg->getStateAlongTrajectory( currentLeg->getLegArrivalTime( ) );
        Eigen::Vector6d nextLegInitialState = nextLeg->getStateAlongTrajectory( nextLeg->getLegDepartureTime() );

        // Add 1 to z position because it's value is approximately zero
        Eigen::Vector3d vectorToAdd (0, 0, 1);
        currentLegFinalState.segment(0, 3) += vectorToAdd;
        nextLegInitialState.segment(0, 3) += vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentLegFinalState.segment(0, 3), nextLegInitialState.segment(0, 3), 1e-5);
    }
}

//! Test delta-V computation for another MGA-1DSM Velocity Formulation trajectory model.
BOOST_AUTO_TEST_CASE( testMGA1DSMVFTrajectory2 )
{
    // Expected test result based on the ideal Cassini 2 trajectory as modelled by GTOP software
    // distributed and downloadable from the ESA website, or within the PaGMO Astrotoolbox.
    const double expectedDeltaV = 8385.15784516116;

    // Create environment
    simulation_setup::SystemOfBodies bodies = createSimplifiedSystemOfBodies( );

    // Set transfer order
    std::vector< std::string > bodyOrder = { "Earth", "Venus", "Venus",  "Earth", "Jupiter", "Saturn" };
    int numberOfNodes = bodyOrder.size( );

    // Create leg settings
    std::vector< std::shared_ptr< TransferLegSettings > > transferLegSettings;
    transferLegSettings.resize( numberOfNodes - 1 );
    transferLegSettings[ 0 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 1 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 2 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 3 ] = dsmVelocityBasedLeg( );
    transferLegSettings[ 4 ] = dsmVelocityBasedLeg( );


    std::vector< std::shared_ptr< TransferNodeSettings > > transferNodeSettings;
    transferNodeSettings.resize( numberOfNodes );
    transferNodeSettings[ 0 ] = escapeAndDepartureNode( std::numeric_limits< double >::infinity( ), 0.0 );
    transferNodeSettings[ 1 ] = swingbyNode( );
    transferNodeSettings[ 2 ] = swingbyNode( );
    transferNodeSettings[ 3 ] = swingbyNode( );
    transferNodeSettings[ 4 ] = swingbyNode( );
    transferNodeSettings[ 5 ] = captureAndInsertionNode( std::numeric_limits< double >::infinity( ), 0.0 );


    std::shared_ptr< TransferTrajectory > transferTrajectory = createTransferTrajectory(
                bodies, transferLegSettings, transferNodeSettings, bodyOrder, "Sun" );

    // Add the time of flight and start epoch, which are in JD.
    double JD = physical_constants::JULIAN_DAY;
    std::vector< double > nodeTimes;
    nodeTimes.push_back( ( -779.046753814506 - 0.5 ) * JD );
    nodeTimes.push_back( nodeTimes.at( 0 ) + 167.378952534645 * JD );
    nodeTimes.push_back( nodeTimes.at( 1 ) + 424.028254165204 * JD );
    nodeTimes.push_back( nodeTimes.at( 2 ) + 53.2897409769205  * JD );
    nodeTimes.push_back( nodeTimes.at( 3 ) + 589.766954923325 * JD );
    nodeTimes.push_back( nodeTimes.at( 4 ) + 2200.00000000000 * JD );

    std::vector< Eigen::VectorXd > transferLegFreeParameters;
    transferLegFreeParameters.resize( numberOfNodes - 1 );
    transferLegFreeParameters[ 0 ] = ( Eigen::VectorXd( 1 ) << 0.769483451363201 ).finished( );
    transferLegFreeParameters[ 1 ] = ( Eigen::VectorXd( 1 ) << 0.513289529822621 ).finished( );
    transferLegFreeParameters[ 2 ] = ( Eigen::VectorXd( 1 ) << 0.0274175362264024 ).finished( );
    transferLegFreeParameters[ 3 ] = ( Eigen::VectorXd( 1 ) << 0.263985256705873 ).finished( );
    transferLegFreeParameters[ 4 ] = ( Eigen::VectorXd( 1 ) << 0.599984695281461 ).finished( );

    std::vector< Eigen::VectorXd > transferNodeFreeParameters;
    transferNodeFreeParameters.resize( numberOfNodes );
    transferNodeFreeParameters[ 0 ] = ( Eigen::VectorXd( 3 ) <<3259.11446832345, 0.525976214695235 * 2 * 3.14159265358979, std::acos(  2 * 0.38086496458657 - 1 ) - 3.14159265358979 / 2 ).finished( );; // 1st leg.
    transferNodeFreeParameters[ 1 ] = ( Eigen::VectorXd( 3 ) <<1.34877968657176 * 6.052e6, -1.5937371121191, 0.0  ).finished( ); // 2nd leg.
    transferNodeFreeParameters[ 2 ] = ( Eigen::VectorXd( 3 ) <<1.05 * 6.052e6, -1.95952512232447, 0.0  ).finished( );// 3rd leg.
    transferNodeFreeParameters[ 3 ] = ( Eigen::VectorXd( 3 ) <<1.30730278372017 * 6.378e6, -1.55498859283059, 0.0  ).finished( ); //4th leg.
    transferNodeFreeParameters[ 4 ] = ( Eigen::VectorXd( 3 ) <<69.8090142993495 * 7.1492e7, -1.5134625299674, 0.0  ).finished( );//4th leg.
    transferNodeFreeParameters[ 5 ] = Eigen::VectorXd( 0 );


    transferTrajectory->evaluateTrajectory(
                nodeTimes, transferLegFreeParameters, transferNodeFreeParameters );
    std::cout << "Transfer with multiple velocity-based DSM legs:" << std::endl;
    printTransferParameterDefinition( transferLegSettings, transferNodeSettings );


    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, transferTrajectory->getTotalDeltaV( ), 1.0E-3);

    // Check continuity of retrieved position between legs
    // Tolerance is large because the state is retrieved from the converged state history (so not from the leg's
    // boundary conditions)
    for( int i = 0; i < transferTrajectory->getNumberOfLegs() - 1; i++ )
    {
        std::shared_ptr< TransferLeg > currentLeg = transferTrajectory->getLegs().at(i);
        std::shared_ptr< TransferLeg > nextLeg = transferTrajectory->getLegs().at(i+1);

        Eigen::Vector6d currentLegFinalState = currentLeg->getStateAlongTrajectory( currentLeg->getLegArrivalTime( ) );
        Eigen::Vector6d nextLegInitialState = nextLeg->getStateAlongTrajectory( nextLeg->getLegDepartureTime() );

        // Add 1 to z position because it's value is approximately zero
        Eigen::Vector3d vectorToAdd (0, 0, 1);
        currentLegFinalState.segment(0, 3) += vectorToAdd;
        nextLegInitialState.segment(0, 3) += vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentLegFinalState.segment(0, 3), nextLegInitialState.segment(0, 3), 1e-5);
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
