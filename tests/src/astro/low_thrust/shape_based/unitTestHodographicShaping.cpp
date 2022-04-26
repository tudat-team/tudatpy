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
 *      Wakker, K. F. (2007), Lecture Notes astro II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h"
#include "tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h"
#include "tudat/astro/low_thrust/shape_based/createBaseFunctionHodographicShaping.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

namespace tudat
{
namespace unit_tests
{

using namespace shape_based_methods;

//! Test hodographic shaping implementation.
BOOST_AUTO_TEST_SUITE( test_hodographic_shaping )

//! Test Earth-Mars transfers, based on Gondelach, D. J., A Hodographic-Shaping Method for Low-Thrust Trajectory Design,
//! TU Delft (MSc thesis), 2012
BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mars_transfer_1 )
{
    spice_interface::loadStandardSpiceKernels( );

    // Basic settings
    int numberOfRevolutions = 2;
    double julianDate = 9264.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1070.0;
    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

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

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
                "Earth"  );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
                "Mars"  );
    Eigen::Vector6d cartesianStateDepartureBody =
            pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Define departure and arrival velocity functions
    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                departureVelocityFunction, arrivalVelocityFunction,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
    hodographicShapingLeg.updateLegParameters(
            ( Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlight * physical_constants::JULIAN_DAY,
            numberOfRevolutions).finished( ) );

    // Compute peak thrust acceleration
    double peakThrustAcceleration = 0.0;
    std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
    hodographicShapingLeg.getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

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
    double expectedDeltaV = 7751.0;
    double expectedPeakAcceleration = 2.64e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  hodographicShapingLeg.getLegDeltaV( ) - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-6 );

}

BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mars_transfer_2 )
{

    // Basic settings
    int numberOfRevolutions = 2;
    double julianDate = 10034.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1070.0;
    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >
            ( ( numberOfRevolutions + 0.5 ) * frequency );
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

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
                "Earth"  );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
                "Mars"  );
    Eigen::Vector6d cartesianStateDepartureBody =
            pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Define departure and arrival velocity functions
    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                departureVelocityFunction, arrivalVelocityFunction,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
    hodographicShapingLeg.updateLegParameters(
            ( Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlight * physical_constants::JULIAN_DAY,
            numberOfRevolutions).finished( ) );

    // Compute peak thrust acceleration
    double peakThrustAcceleration = 0.0;
    std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
    hodographicShapingLeg.getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

    for( auto it : thrustAccelerationsAlongTrajectory )
    {
        Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
        if ( currentCartesianThrustAcceleration.norm() > peakThrustAcceleration )
        {
            peakThrustAcceleration = currentCartesianThrustAcceleration.norm();
        }
    }

    // Check results consistency w.r.t. Gondelach, D. J., A Hodographic-Shaping Method for Low-Thrust Trajectory Design,
    //! TU Delft (MSc thesis), 2012
    double expectedDeltaV = 6742.0;
    double expectedPeakAcceleration = 1.46e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  hodographicShapingLeg.getLegDeltaV( ) - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-6 );
}

BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mars_transfer_3 )
{

    // Basic settings
    int numberOfRevolutions = 2;
    double julianDate = 9244.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1090.0;
    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

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

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
                "Earth"  );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
                "Mars"  );
    Eigen::Vector6d cartesianStateDepartureBody =
            pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Define departure and arrival velocity functions
    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                departureVelocityFunction, arrivalVelocityFunction,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
    hodographicShapingLeg.updateLegParameters(
            ( Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlight * physical_constants::JULIAN_DAY,
            numberOfRevolutions).finished( ) );

    // Compute peak thrust acceleration
    double peakThrustAcceleration = 0.0;
    std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
    hodographicShapingLeg.getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

    for( auto it : thrustAccelerationsAlongTrajectory )
    {
        Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
        if ( currentCartesianThrustAcceleration.norm() > peakThrustAcceleration )
        {
            peakThrustAcceleration = currentCartesianThrustAcceleration.norm();
        }
    }

    // Check results consistency w.r.t. Gondelach, D. J., A Hodographic-Shaping Method for Low-Thrust Trajectory Design,
    //! TU Delft (MSc thesis), 2012
    double expectedDeltaV = 6686.0;
    double expectedPeakAcceleration = 2.46e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  hodographicShapingLeg.getLegDeltaV( ) - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-6 );
}

BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mars_transfer_4 )
{

    // Basic settings
    int numberOfRevolutions = 2;
    double julianDate = 10024.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1050.0;
    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );

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

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
                "Earth"  );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
                "Mars"  );
    Eigen::Vector6d cartesianStateDepartureBody =
            pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Define departure and arrival velocity functions
    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                departureVelocityFunction, arrivalVelocityFunction,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
    hodographicShapingLeg.updateLegParameters(
            ( Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlight * physical_constants::JULIAN_DAY,
            numberOfRevolutions).finished( ) );

    // Compute peak thrust acceleration
    double peakThrustAcceleration = 0.0;
    std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
    hodographicShapingLeg.getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

    for( auto it : thrustAccelerationsAlongTrajectory )
    {
        Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
        if ( currentCartesianThrustAcceleration.norm() > peakThrustAcceleration )
        {
            peakThrustAcceleration = currentCartesianThrustAcceleration.norm();
        }
    }

    // Check results consistency w.r.t. Gondelach, D. J., A Hodographic-Shaping Method for Low-Thrust Trajectory Design,
    //! TU Delft (MSc thesis), 2012
    double expectedDeltaV = 6500.0;
    double expectedPeakAcceleration = 1.58e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  hodographicShapingLeg.getLegDeltaV( ) - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-6 );
}

BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mars_transfer_5 )
{
    // Basic settings
    int numberOfRevolutions = 2;
    double julianDate = 10024.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 1050.0;
    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the normal velocity composite function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );

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

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
                "Earth"  );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
                "Mars"  );
    Eigen::Vector6d cartesianStateDepartureBody =
            pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Define departure and arrival velocity functions
    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                departureVelocityFunction, arrivalVelocityFunction,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
    hodographicShapingLeg.updateLegParameters(
            ( Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlight * physical_constants::JULIAN_DAY,
            numberOfRevolutions).finished( ) );

    // Compute peak thrust acceleration
    double peakThrustAcceleration = 0.0;
    std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
    hodographicShapingLeg.getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

    for( auto it : thrustAccelerationsAlongTrajectory )
    {
        Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
        if ( currentCartesianThrustAcceleration.norm() > peakThrustAcceleration )
        {
            peakThrustAcceleration = currentCartesianThrustAcceleration.norm();
        }
    }

    // Check results consistency w.r.t. Gondelach, D. J., A Hodographic-Shaping Method for Low-Thrust Trajectory Design,
    //! TU Delft (MSc thesis), 2012
    double expectedDeltaV = 6342.0;
    double expectedPeakAcceleration = 1.51e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  hodographicShapingLeg.getLegDeltaV( ) - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-6 );

}


//! Test Earth-Mercury transfer, based on the thesis by Gondelach.
BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mercury_transfer )
{
    using namespace shape_based_methods;

    /// First Earth-Mercury transfer.

    int numberOfRevolutions = 1;
    double julianDate = 5025 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 440.0;
    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

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

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 6.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                6.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
                "Earth"  );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
                "Mercury" );
    Eigen::Vector6d cartesianStateDepartureBody =
            pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Define departure and arrival velocity functions
    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                departureVelocityFunction, arrivalVelocityFunction,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
    hodographicShapingLeg.updateLegParameters(
            ( Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlight * physical_constants::JULIAN_DAY,
            numberOfRevolutions).finished( ) );

    // Compute peak thrust acceleration
    double peakThrustAcceleration = 0.0;
    std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
    hodographicShapingLeg.getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

    for( auto it : thrustAccelerationsAlongTrajectory )
    {
        Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
        if ( currentCartesianThrustAcceleration.norm() > peakThrustAcceleration )
        {
            peakThrustAcceleration = currentCartesianThrustAcceleration.norm();
        }
    }

    // Check results consistency w.r.t. Gondelach, D. J., A Hodographic-Shaping Method for Low-Thrust Trajectory Design,
    //! TU Delft (MSc thesis), 2012
    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in the used Mercury ephemeris.
    double expectedDeltaV = 28082.0;
    double expectedPeakAcceleration = 64.1e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  hodographicShapingLeg.getLegDeltaV( ) - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-6 );
}

//! Test Earth-Mars transfer with shaping functions with free coefficients. Compare with the DeltaV calculated using the
//! original version of the hodographic shaping code.
BOOST_AUTO_TEST_CASE( test_hodographic_shaping_with_free_parameters )
{

    double numberOfRevolutions = 1.0;
    double julianDate = 2458849.5; // julian date specified in SECONDS
    double timeOfFlight = 500.0; // time of flight specified in DAYS
    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

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

    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = ( Eigen::Vector2d( ) << 500.0, 500.0 ).finished( );
    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = ( Eigen::Vector2d( ) << 500.0, -200.0 ).finished( );
    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = ( Eigen::Vector2d( ) << 500.0, 2000.0 ).finished( );

    // Retrieve cartesian state at departure and arrival.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
                "Earth"  );
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
                "Mars"  );
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * physical_constants::JULIAN_DAY );

    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };

    // Define departure and arrival velocity functions
    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                departureVelocityFunction, arrivalVelocityFunction,
                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
    hodographicShapingLeg.updateLegParameters(
            ( Eigen::VectorXd( 9 ) << julianDate, julianDate + timeOfFlight  * physical_constants::JULIAN_DAY,
            numberOfRevolutions,  freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
            freeCoefficientsAxialVelocityFunction ).finished( ) );

    // Compare computed deltaV with deltaV predicted by original version of Gondelach's code.
    double expectedDeltaV = 172313.0;
    BOOST_CHECK_SMALL(std::fabs(hodographicShapingLeg.getLegDeltaV( ) - expectedDeltaV), 1.0);
}

//    /// Second Earth-Mercury transfer.

//    numberOfRevolutions = 1;

//    julianDate = 5015 * physical_constants::JULIAN_DAY;

//    timeOfFlight = 450.0;

//    // Set vehicle mass.
//    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

//    // Define integrator settings.
//    integratorSettings = std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                numerical_integrators::rungeKutta4, 0.0, timeOfFlight * physical_constants::JULIAN_DAY / 500.0 );

//    // Retrieve cartesian state at departure and arrival.
//    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
//    cartesianStateArrivalBody =
//            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

//    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    // Create base function settings for the components of the radial velocity composite function.
//    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

//    // Create components of the radial velocity composite function.
//    radialVelocityFunctionComponents.clear( );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

//    // Create base function settings for the components of the normal velocity composite function.
//    firstNormalVelocityBaseFunctionSettings =
//            std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdNormalVelocityBaseFunctionSettings =
//            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

//    // Create components of the normal velocity composite function.
//    normalVelocityFunctionComponents.clear( );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );


//    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//    hodographicShaping = HodographicShaping(
//                cartesianStateDepartureBody, cartesianStateArrivalBody,
//                timeOfFlight * physical_constants::JULIAN_DAY, numberOfRevolutions,
//                bodies, "Vehicle", "Sun",
//                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
//                integratorSettings );

//    stepSize = ( timeOfFlight * physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
//    peakAcceleration = 0.0;

//    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
//        double currentTime = currentStep * stepSize;

//        double currentAcceleration = hodographicShaping.computeThrustAccelerationVector( currentTime ).norm( );
//        if ( currentAcceleration > peakAcceleration )
//        {
//            peakAcceleration = currentAcceleration;
//        }

//    }


//    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
//    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
//    expectedDeltaV = 26997.0;
//    expectedPeakAcceleration = 63.5e-4;

//    // DeltaV provided with a precision of 1 m/s
//    BOOST_CHECK_SMALL( std::fabs(  hodographicShaping.computeDeltaV( ) - expectedDeltaV ), 100.0 );
//    // Peak acceleration provided with a precision 1.0e-6 m/s^2
//    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );


//    /// Third Earth-Mercury transfer.

//    numberOfRevolutions = 0;

//    julianDate = 4675 * physical_constants::JULIAN_DAY;

//    timeOfFlight = 190.0;

//    // Set vehicle mass.
//    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

//    // Define integrator settings.
//    integratorSettings = std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                numerical_integrators::rungeKutta4, 0.0, timeOfFlight * physical_constants::JULIAN_DAY / 500.0 );

//    // Retrieve cartesian state at departure and arrival.
//    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
//    cartesianStateArrivalBody =
//            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

//    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    // Create base function settings for the components of the radial velocity composite function.
//    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

//    // Create components of the radial velocity composite function.
//    radialVelocityFunctionComponents.clear( );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

//    // Create base function settings for the components of the normal velocity composite function.
//    firstNormalVelocityBaseFunctionSettings =
//            std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdNormalVelocityBaseFunctionSettings =
//            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency * 0.5 );

//    // Create components of the normal velocity composite function.
//    normalVelocityFunctionComponents.clear( );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );


//    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//    hodographicShaping = HodographicShaping(
//                cartesianStateDepartureBody, cartesianStateArrivalBody,
//                timeOfFlight * physical_constants::JULIAN_DAY, numberOfRevolutions,
//                bodies, "Vehicle", "Sun",
//                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
//                integratorSettings );

//    stepSize = ( timeOfFlight * physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
//    peakAcceleration = 0.0;

//    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
//        double currentTime = currentStep * stepSize;

//        double currentAcceleration = hodographicShaping.computeThrustAccelerationVector( currentTime ).norm( );
//        if ( currentAcceleration > peakAcceleration )
//        {
//            peakAcceleration = currentAcceleration;
//        }

//    }


//    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
//    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
//    expectedDeltaV = 22683.0;
//    expectedPeakAcceleration = 56.0e-4;

//    // DeltaV provided with a precision of 1 m/s
//    BOOST_CHECK_SMALL( std::fabs(  hodographicShaping.computeDeltaV( ) - expectedDeltaV ), 100.0 );
//    // Peak acceleration provided with a precision 1.0e-6 m/s^2
//    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );



//    /// Fourth Earth-Mercury transfer.

//    numberOfRevolutions = 0;

//    julianDate = 4675 * physical_constants::JULIAN_DAY;

//    timeOfFlight = 190.0;

//    // Set vehicle mass.
//    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

//    // Define integrator settings.
//    integratorSettings = std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                numerical_integrators::rungeKutta4, 0.0, timeOfFlight * physical_constants::JULIAN_DAY / 500.0 );

//    // Retrieve cartesian state at departure and arrival.
//    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
//    cartesianStateArrivalBody =
//            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

//    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    // Create base function settings for the components of the radial velocity composite function.
//    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdRadialVelocityBaseFunctionSettings = std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency * 0.5 );

//    // Create components of the radial velocity composite function.
//    radialVelocityFunctionComponents.clear( );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( sine, thirdRadialVelocityBaseFunctionSettings ) );

//    // Create base function settings for the components of the normal velocity composite function.
//    firstNormalVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondNormalVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdNormalVelocityBaseFunctionSettings = std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

//    // Create components of the normal velocity composite function.
//    normalVelocityFunctionComponents.clear( );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );


//    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//    hodographicShaping = HodographicShaping(
//                cartesianStateDepartureBody, cartesianStateArrivalBody,
//                timeOfFlight * physical_constants::JULIAN_DAY, numberOfRevolutions,
//                bodies, "Vehicle", "Sun",
//                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
//                integratorSettings );

//    stepSize = ( timeOfFlight * physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
//    peakAcceleration = 0.0;

//    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
//        double currentTime = currentStep * stepSize;

//        double currentAcceleration = hodographicShaping.computeThrustAccelerationVector( currentTime ).norm( );
//        if ( currentAcceleration > peakAcceleration )
//        {
//            peakAcceleration = currentAcceleration;
//        }

//    }


//    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
//    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
//    expectedDeltaV = 22613.0;
//    expectedPeakAcceleration = 56.0e-4;

//    // DeltaV provided with a precision of 1 m/s
//    BOOST_CHECK_SMALL( std::fabs(  hodographicShaping.computeDeltaV( ) - expectedDeltaV ), 100.0 );
//    // Peak acceleration provided with a precision 1.0e-6 m/s^2
//    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );



//    /// Fifth Earth-Mercury transfer.

//    numberOfRevolutions = 0;

//    julianDate = 4355 * physical_constants::JULIAN_DAY;

//    timeOfFlight = 160.0;

//    // Set vehicle mass.
//    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

//    // Define integrator settings.
//    integratorSettings = std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                numerical_integrators::rungeKutta4, 0.0, timeOfFlight * physical_constants::JULIAN_DAY / 500.0 );

//    // Retrieve cartesian state at departure and arrival.
//    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
//    cartesianStateArrivalBody =
//            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

//    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    // Create base function settings for the components of the radial velocity composite function.
//    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

//    // Create components of the radial velocity composite function.
//    radialVelocityFunctionComponents.clear( );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

//    // Create base function settings for the components of the normal velocity composite function.
//    firstNormalVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    secondNormalVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    thirdNormalVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

//    // Create components of the normal velocity composite function.
//    normalVelocityFunctionComponents.clear( );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );


//    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//    hodographicShaping = HodographicShaping(
//                cartesianStateDepartureBody, cartesianStateArrivalBody,
//                timeOfFlight * physical_constants::JULIAN_DAY, numberOfRevolutions,
//                bodies, "Vehicle", "Sun",
//                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
//                integratorSettings );

//    stepSize = ( timeOfFlight * physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
//    peakAcceleration = 0.0;

//    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
//        double currentTime = currentStep * stepSize;

//        double currentAcceleration = hodographicShaping.computeThrustAccelerationVector( currentTime ).norm( );
//        if ( currentAcceleration > peakAcceleration )
//        {
//            peakAcceleration = currentAcceleration;
//        }

//    }


//    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
//    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
//    expectedDeltaV = 21766.0;
//    expectedPeakAcceleration = 50.3e-4;

//    // DeltaV provided with a precision of 1 m/s
//    BOOST_CHECK_SMALL( std::fabs(  hodographicShaping.computeDeltaV( ) - expectedDeltaV ), 100.0 );
//    // Peak acceleration provided with a precision 1.0e-6 m/s^2
//    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );

//}

//SystemOfBodies getTestBodyMap( )
//{
//    spice_interface::loadStandardSpiceKernels( );
//
//    // Create central, departure and arrival bodies.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Sun" );
//    bodiesToCreate.push_back( "Earth" );
//    bodiesToCreate.push_back( "Mars" );
//    bodiesToCreate.push_back( "Jupiter" );
//
//
//    std::string frameOrigin = "SSB";
//    std::string frameOrientation = "ECLIPJ2000";
//
//    BodyListSettings bodySettings =
//            getDefaultBodySettings( bodiesToCreate, frameOrigin, frameOrientation );
//
//
//
//    // Define central body ephemeris settings.
//    bodySettings.at( "Sun" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
//                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );
//
//
//    // Create system of bodies.
//    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
//
//    bodies.createEmptyBody( "Vehicle" );
//    return bodies;
//}
//
////! Test.
//BOOST_AUTO_TEST_CASE( test_hodographic_shaping_full_propagation )
//{
//
//    double numberOfRevolutions = 1.0;
//    double julianDate = 2458849.5; // julian date specified in SECONDS
//    double timeOfFlight = 500.0; // time of flight specified in DAYS
//    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
//    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );
//
//    // Create base function settings for the components of the radial velocity composite function.
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
//            std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
//
//    // Create components of the radial velocity composite function.
//    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );
//
//    // Create base function settings for the components of the normal velocity composite function.
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
//            std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
//
//    // Create components of the normal velocity composite function.
//    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, fourthNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthNormalVelocityBaseFunctionSettings ) );
//
//    // Create base function settings for the components of the axial velocity composite function.
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
//            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
//            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
//
//    // Set components for the axial velocity function.
//    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, fourthAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, fifthAxialVelocityBaseFunctionSettings ) );
//
//    // Initialize free coefficients vector for radial velocity function.
//    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = ( Eigen::Vector2d( ) << 500.0, 500.0 ).finished( );
//    // Initialize free coefficients vector for normal velocity function.
//    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = ( Eigen::Vector2d( ) << 500.0, -200.0 ).finished( );
//    // Initialize free coefficients vector for axial velocity function.
//    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = ( Eigen::Vector2d( ) << 500.0, 2000.0 ).finished( );
//
//    // Retrieve cartesian state at departure and arrival.
//    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
//                "Earth"  );
//    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
//                "Mars"  );
//    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
//    Eigen::Vector6d cartesianStateArrivalBody =
//            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * physical_constants::JULIAN_DAY );
//
//    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//    std::shared_ptr< HodographicShaping > hodographicShaping = std::make_shared< HodographicShaping >(
//                cartesianStateDepartureBody, cartesianStateArrivalBody,
//                timeOfFlight * physical_constants::JULIAN_DAY,
//                spice_interface::getBodyGravitationalParameter( "Sun" ), numberOfRevolutions,
//                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction );
//
//    std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return cartesianStateDepartureBody.segment( 3, 3 ); };
//    std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return cartesianStateArrivalBody.segment( 3, 3 ); };
//
//    // Define departure and arrival velocity functions
//    HodographicShapingLeg hodographicShapingLeg = HodographicShapingLeg(
//                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
//                spice_interface::getBodyGravitationalParameter( "Sun" ),
//                departureVelocityFunction, arrivalVelocityFunction,
//                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents );
//    hodographicShapingLeg.updateLegParameters(
//            ( Eigen::VectorXd( 9 ) << julianDate, julianDate + timeOfFlight  * physical_constants::JULIAN_DAY,
//            numberOfRevolutions,  freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
//            freeCoefficientsAxialVelocityFunction ).finished( ) );
//
//    // Create environment
//    SystemOfBodies bodies = getTestBodyMap( );
//
//    // Define mass function of the vehicle.
//    std::function< double( const double ) > newMassFunction = [ = ]( const double ){ return 2000.0; }; // - 50.0 / ( timeOfFlight * physical_constants::JULIAN_DAY ) * currentTime ;
//    bodies.at( "Vehicle" )->setBodyMassFunction( newMassFunction );
//
//    // Define integrator settings.
//    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
//            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
//                numerical_integrators::rungeKutta4, 0.0,  timeOfFlight * physical_constants::JULIAN_DAY / 1000.0 );
//
//    // Create object with list of dependent variables
//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
//    dependentVariablesList.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
//                                          basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
//    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
//            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );
//
//    // Create termination conditions settings.
//    std::pair< std::shared_ptr< PropagationTerminationSettings >,
//            std::shared_ptr< PropagationTerminationSettings > > terminationConditions = std::make_pair(
//                std::make_shared< PropagationTimeTerminationSettings >( 0.0 ),
//                std::make_shared< PropagationTimeTerminationSettings >( timeOfFlight * physical_constants::JULIAN_DAY ) );
//
//    // Create complete propagation settings (backward and forward propagations).
//    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double ){ return 3000.0; };
//    basic_astrodynamics::AccelerationMap lowThrustAccelerationsMap = retrieveLowThrustAccelerationMap(
//                   hodographicShaping, bodies, "Vehicle", "Sun", specificImpulseFunction, 0.0 );
//
//
//    std::pair< std::shared_ptr< PropagatorSettings< double > >,
//            std::shared_ptr< PropagatorSettings< double > > > propagatorSettings =
//            createLowThrustTranslationalStatePropagatorSettings(
//                hodographicShaping, "Vehicle", "Sun", lowThrustAccelerationsMap, dependentVariablesToSave );
//
//
//    // Compute shaped trajectory and propagated trajectory.
//    std::map< double, Eigen::VectorXd > fullPropagationResults;
//    std::map< double, Eigen::Vector6d > shapingMethodResults;
//    std::map< double, Eigen::VectorXd > dependentVariablesHistory;
//
//    computeLowThrustLegSemiAnalyticalAndFullPropagation(
//                hodographicShaping, bodies, integratorSettings, propagatorSettings,
//                fullPropagationResults, shapingMethodResults, dependentVariablesHistory );
//
//    // Check that boundary conditions are still fulfilled when free parameters are added.
//    for ( int i = 0 ; i < 6 ; i++ )
//    {
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin( )->second[ i ] - cartesianStateDepartureBody[ i ] )
//                           / shapingMethodResults.begin( )->second[ i ], 1.0e-8 );
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin( )->second[ i ] - cartesianStateArrivalBody[ i ] )
//                           / shapingMethodResults.rbegin( )->second[ i ], 1.0e-8 );
//    }
//
//    // Check results consistency between full propagation and shaped trajectory at departure and arrival.
//    for ( int i = 0 ; i < 6 ; i++ )
//    {
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin( )->second[ i ] - fullPropagationResults.begin( )->second[ i ] )
//                           / shapingMethodResults.begin( )->second[ i ], 2.0e-7 );
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin( )->second[ i ] - fullPropagationResults.rbegin( )->second[ i ] )
//                           / shapingMethodResults.rbegin( )->second[ i ], 2.0e-7 );
//    }
//
//}


////! Test full propagation while propagating the spacecraft mass too.
//BOOST_AUTO_TEST_CASE( test_hodographic_shaping_full_propagation_mass_propagation )
//{
//    double numberOfRevolutions = 1.0;
//    double julianDate = 2458849.5;
//    double timeOfFlight = 500.0;
//    double initialBodyMass = 2000.0;
//    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * physical_constants::JULIAN_DAY );
//    double scaleFactor = 1.0 / ( timeOfFlight * physical_constants::JULIAN_DAY );

//    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double )
//    { return 3000.0; };

//    // Retrieve cartesian state at departure and arrival.
//    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris>(
//                "Earth"  );
//    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximateJplEphemeris >(
//                "Mars"  );
//    Eigen::Vector6d cartesianStateDepartureBody =
//            pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
//    Eigen::Vector6d cartesianStateArrivalBody =
//            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * physical_constants::JULIAN_DAY );

//    // Create base function settings for the components of the radial velocity composite function.
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
//            std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthRadialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );

//    // Create components of the radial velocity composite function.
//    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > radialVelocityFunctionComponents;
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, fourthRadialVelocityBaseFunctionSettings ) );
//    radialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthRadialVelocityBaseFunctionSettings ) );

//    // Create base function settings for the components of the normal velocity composite function.
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
//            std::make_shared< BaseFunctionHodographicShapingSettings >( );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthNormalVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >( 1.0, 0.5 * frequency, scaleFactor );

//    // Create components of the normal velocity composite function.
//    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, fourthNormalVelocityBaseFunctionSettings ) );
//    normalVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, fifthNormalVelocityBaseFunctionSettings ) );

//    // Create base function settings for the components of the axial velocity composite function.
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
//            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
//            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fourthAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
//    std::shared_ptr< BaseFunctionHodographicShapingSettings > fifthAxialVelocityBaseFunctionSettings =
//            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
//                4.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

//    // Set components for the axial velocity function.
//    std::vector< std::shared_ptr< BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerCosine, fourthAxialVelocityBaseFunctionSettings ) );
//    axialVelocityFunctionComponents.push_back(
//                createBaseFunctionHodographicShaping( scaledPowerSine, fifthAxialVelocityBaseFunctionSettings ) );

//    // Initialize free coefficients vector for radial velocity function.
//    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
//    freeCoefficientsRadialVelocityFunction[ 0 ] = 500.0;
//    freeCoefficientsRadialVelocityFunction[ 1 ] = 500.0;

//    // Initialize free coefficients vector for normal velocity function.
//    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 2 );
//    freeCoefficientsNormalVelocityFunction[ 0 ] = 500.0;
//    freeCoefficientsNormalVelocityFunction[ 1 ] = - 200.0;

//    // Initialize free coefficients vector for axial velocity function.
//    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 2 );
//    freeCoefficientsAxialVelocityFunction[ 0 ] = 500.0;
//    freeCoefficientsAxialVelocityFunction[ 1 ] = 2000.0;


//    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
//    HodographicShaping hodographicShaping(
//                cartesianStateDepartureBody, cartesianStateArrivalBody,
//                timeOfFlight * physical_constants::JULIAN_DAY,
//                spice_interface::getBodyGravitationalParameter( "Sun" ), 1,
//                radialVelocityFunctionComponents, normalVelocityFunctionComponents, axialVelocityFunctionComponents,
//                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction, freeCoefficientsAxialVelocityFunction,
//                initialBodyMass );

//    // Create environment
//    SystemOfBodies bodies = getTestBodyMap( );
//    bodies.at( "Vehicle" )->setConstantBodyMass( initialBodyMass );

//    // Define integrator settings.
//    double stepSize = ( timeOfFlight * physical_constants::JULIAN_DAY ) / static_cast< double >( 50 );
//    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
//            std::make_shared< numerical_integrators::IntegratorSettings< double > >(
//                numerical_integrators::rungeKutta4, 0.0, stepSize / 400.0 );


//    // Define list of dependent variables to save.
//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
//    dependentVariablesList.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
//                                          basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
//    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
//                                          total_mass_rate_dependent_variables, "Vehicle" ) );

//    // Create object with list of dependent variables
//    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
//            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );

//    // Create termination conditions settings.
//    std::pair< std::shared_ptr< PropagationTerminationSettings >, std::shared_ptr< PropagationTerminationSettings > > terminationConditions;
//    terminationConditions.first = std::make_shared< PropagationTimeTerminationSettings >( 0.0 );
//    terminationConditions.second = std::make_shared< PropagationTimeTerminationSettings >( timeOfFlight * physical_constants::JULIAN_DAY );

//    // Create complete propagation settings (backward and forward propagations).
//    std::pair< std::shared_ptr< PropagatorSettings< double > >,
//            std::shared_ptr< PropagatorSettings< double > > > propagatorSettings = hodographicShaping.createLowThrustPropagatorSettings(
//                bodies, "Vehicle", "Sun", specificImpulseFunction, basic_astrodynamics::AccelerationMap( ), integratorSettings,
//                dependentVariablesToSave );

//    // Compute shaped trajectory and propagated trajectory.
//    std::map< double, Eigen::VectorXd > fullPropagationResults;
//    std::map< double, Eigen::Vector6d > shapingMethodResults;
//    std::map< double, Eigen::VectorXd > dependentVariablesHistory;
//    hodographicShaping.computeSemiAnalyticalAndFullPropagation(
//                bodies, integratorSettings, propagatorSettings,
//                fullPropagationResults, shapingMethodResults, dependentVariablesHistory );

//    // Check that boundary conditions are still fulfilled when free parameters are added.
//    for ( int i = 0 ; i < 6 ; i++ )
//    {
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin( )->second[ i ] - cartesianStateDepartureBody[ i ] )
//                           / shapingMethodResults.begin( )->second[ i ], 1.0e-8 );
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin( )->second[ i ] - cartesianStateArrivalBody[ i ] )
//                           / shapingMethodResults.rbegin( )->second[ i ], 1.0e-8 );
//    }

//    // Check results consistency between full propagation and shaped trajectory at departure and arrival.
//    for ( int i = 0 ; i < 6 ; i++ )
//    {
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin( )->second[ i ] - fullPropagationResults.begin( )->second[ i ] )
//                           / shapingMethodResults.begin( )->second[ i ], 2.0e-7 );
//        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin( )->second[ i ] - fullPropagationResults.rbegin( )->second[ i ] )
//                           / shapingMethodResults.rbegin( )->second[ i ], 2.0e-7 );
//    }

//    // Check consistency between current and expected mass rates.
//    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin( ) ; itr != dependentVariablesHistory.end( ) ; itr++ )
//    {
//        Eigen::Vector3d currentThrustVector = itr->second.segment( 0, 3 );
//        double currentMass = fullPropagationResults.at( itr->first )( 6 );
//        double currentMassRate = - itr->second( 3 );
//        double expectedMassRate = currentThrustVector.norm( ) * currentMass /
//                ( specificImpulseFunction( itr->first ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

//        BOOST_CHECK_SMALL( std::fabs( currentMassRate - expectedMassRate ), 1.0e-15 );

//    }


//    // Test trajectory function.
//    std::vector< double > epochsVector;
//    epochsVector.push_back( 0.0 );
//    epochsVector.push_back( timeOfFlight / 4.0 * physical_constants::JULIAN_DAY );
//    epochsVector.push_back( timeOfFlight / 2.0 * physical_constants::JULIAN_DAY );
//    epochsVector.push_back( 3.0 * timeOfFlight / 4.0 * physical_constants::JULIAN_DAY );
//    epochsVector.push_back( timeOfFlight * physical_constants::JULIAN_DAY );

//    std::map< double, Eigen::Vector6d > trajectory;
//    std::map< double, Eigen::VectorXd > massProfile;
//    std::map< double, Eigen::VectorXd > thrustProfile;
//    std::map< double, Eigen::VectorXd > thrustAccelerationProfile;

////    hodographicShaping.getTrajectory( epochsVector, trajectory );
////    hodographicShaping.getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );
//////    hodographicShaping.getThrustForceProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );
//////    hodographicShaping.getThrustAccelerationProfile( epochsVector, thrustAccelerationProfile, specificImpulseFunction, integratorSettings );

////    for ( int i = 0 ; i < 3 ; i ++ )
////    {
////        BOOST_CHECK_SMALL( std::fabs( trajectory.begin( )->second[ i ] - cartesianStateDepartureBody[ i ] ), 1.0e-3 );
////        BOOST_CHECK_SMALL( std::fabs( trajectory.begin( )->second[ i + 3 ] - cartesianStateDepartureBody[ i + 3 ] ), 1.0e-10 );
////        BOOST_CHECK_SMALL( std::fabs( trajectory.rbegin( )->second[ i ] - cartesianStateArrivalBody[ i ] ), 1.0e-3 );
////        BOOST_CHECK_SMALL( std::fabs( trajectory.rbegin( )->second[ i + 3 ] - cartesianStateArrivalBody[ i + 3 ] ), 1.0e-10 );
////    }

////    for ( std::map< double, Eigen::Vector6d >::iterator itr = trajectory.begin( ) ; itr != trajectory.end( ) ; itr++ )
////    {
////        Eigen::Vector6d stateVector = hodographicShaping.computeCurrentCartesianState( itr->first );
////        Eigen::Vector3d thrustAccelerationVector = hodographicShaping.computeThrustAcceleration(
////                    itr->first, specificImpulseFunction, integratorSettings );
////        Eigen::Vector3d thrustVector = hodographicShaping.computeCurrentThrustForce(
////                    itr->first, specificImpulseFunction, integratorSettings );
////        double mass = hodographicShaping.computeCurrentMass( itr->first, specificImpulseFunction, integratorSettings );

////        for ( int i = 0 ; i < 3 ; i++ )
////        {
////            BOOST_CHECK_SMALL( std::fabs( itr->second[ i ] - stateVector[ i ] ), 1.0e-6 );
////            BOOST_CHECK_SMALL( std::fabs( itr->second[ i + 3 ] - stateVector[ i + 3 ] ), 1.0e-12 );
////            BOOST_CHECK_SMALL( std::fabs( thrustAccelerationProfile[ itr->first ][ i ] - thrustAccelerationVector[ i ] ), 1.0e-6 );
////            BOOST_CHECK_SMALL( std::fabs( thrustProfile[ itr->first ][ i ] - thrustVector[ i ] ), 1.0e-12 );
////        }
////        BOOST_CHECK_SMALL( std::fabs( massProfile[ itr->first ][ 0 ] - mass ), 1.0e-12 );
////    }
////}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
