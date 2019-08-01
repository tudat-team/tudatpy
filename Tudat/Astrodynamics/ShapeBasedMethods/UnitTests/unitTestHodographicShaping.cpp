/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes Astrodynamics II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <math.h>
#include <iostream>

#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"

namespace tudat
{
namespace unit_tests
{

//! Test hodographic shaping implementation.
BOOST_AUTO_TEST_SUITE( test_hodographic_shaping )

//! Test Earth-Mars transfers, based on the thesis by Gondelach (ADD PROPER REFERENCE).
BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mars_transfer )
{
    using namespace shape_based_methods;

    int numberOfRevolutions = 2;

    /// First Earth-Mars transfer.

    double julianDate = 9264.5 * physical_constants::JULIAN_DAY;

    double timeOfFlight = 1070.0;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );


    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    double scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping VelocityShapingMethod(
                cartesianStateDepartureBody, cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY, numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );


    int numberOfSteps = 500;
    double stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( numberOfSteps );
    double peakAcceleration = 0.0;

    for ( int currentStep = 0 ; currentStep <= numberOfSteps ; currentStep++ ){

        double currentTime = currentStep * stepSize;

        double currentAccelerationMagnitude = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();

        if ( currentAccelerationMagnitude > peakAcceleration )
        {
            peakAcceleration = currentAccelerationMagnitude;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)

    double expectedDeltaV = 7751.0;
    double expectedPeakAcceleration = 2.64e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-6 );



    /// Second Earth-Mars transfer.

    julianDate = 10034.5 * physical_constants::JULIAN_DAY;

    timeOfFlight = 1070.0;

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings = std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
   firstAxialVelocityBaseFunctionSettings = std::make_shared< TrigonometricFunctionHodographicShapingSettings >
           ( ( numberOfRevolutions + 0.5 ) * frequency );
    secondAxialVelocityBaseFunctionSettings = std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    thirdAxialVelocityBaseFunctionSettings = std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    axialVelocityFunctionComponents.clear();
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY, numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    peakAcceleration = 0.0;

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }
    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)

    expectedDeltaV = 6742.0;
    expectedPeakAcceleration = 1.46e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-6 );



    /// Third Earth-Mars transfer.

    julianDate = 9244.5 * physical_constants::JULIAN_DAY;

    timeOfFlight = 1090.0;

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    axialVelocityFunctionComponents.clear();
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    peakAcceleration = 0.0;

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)

    expectedDeltaV = 6686.0;
    expectedPeakAcceleration = 2.46e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-6 );


    /// Fourth Earth-Mars transfer.

    julianDate = 10024.5 * physical_constants::JULIAN_DAY;

    timeOfFlight = 1050.0;

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    axialVelocityFunctionComponents.clear();
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );


    peakAcceleration = 0.0;

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)

    expectedDeltaV = 6500.0;
    expectedPeakAcceleration = 1.58e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-6 );



    /// Fifth Earth-Mars transfer.

    julianDate = 10024.5 * physical_constants::JULIAN_DAY;

    timeOfFlight = 1050.0;

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the axial velocity composite function.
    firstAxialVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    secondAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    axialVelocityFunctionComponents.clear();
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    peakAcceleration = 0.0;

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)

    expectedDeltaV = 6342.0;
    expectedPeakAcceleration = 1.51e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 1.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-6 );


}


//! Test Earth-Mercury transfer, based on the thesis by Gondelach.
BOOST_AUTO_TEST_CASE( test_hodographic_shaping_earth_mercury_transfer )
{
    using namespace shape_based_methods;

    spice_interface::loadStandardSpiceKernels();

    /// First Earth-Mercury transfer.

    int numberOfRevolutions = 1;

    double julianDate = 5025 * physical_constants::JULIAN_DAY;

    double timeOfFlight = 440.0;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mercury );

    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );


    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    double scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    // Initialize free coefficients vector for radial velocity function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping VelocityShapingMethod(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    double stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    double peakAcceleration = 0.0;

    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }
    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
    double expectedDeltaV = 28082.0;
    double expectedPeakAcceleration = 64.1e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 100.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );



    /// Second Earth-Mercury transfer.

    numberOfRevolutions = 1;

    julianDate = 5015 * physical_constants::JULIAN_DAY;

    timeOfFlight = 450.0;

    // Retrieve cartesian state at departure and arrival.
    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( cosine, thirdNormalVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    peakAcceleration = 0.0;

    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
    expectedDeltaV = 26997.0;
    expectedPeakAcceleration = 63.5e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 100.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );


    /// Third Earth-Mercury transfer.

    numberOfRevolutions = 0;

    julianDate = 4675 * physical_constants::JULIAN_DAY;

    timeOfFlight = 190.0;

    // Retrieve cartesian state at departure and arrival.
    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings =
            std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings =
            std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency * 0.5 );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    peakAcceleration = 0.0;

    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
    expectedDeltaV = 22683.0;
    expectedPeakAcceleration = 56.0e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 100.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );



    /// Fourth Earth-Mercury transfer.

    numberOfRevolutions = 0;

    julianDate = 4675 * physical_constants::JULIAN_DAY;

    timeOfFlight = 190.0;

    // Retrieve cartesian state at departure and arrival.
    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< TrigonometricFunctionHodographicShapingSettings >( frequency * 0.5 );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings = std::make_shared< TrigonometricFunctionHodographicShapingSettings >( 0.5 * frequency );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( sine, thirdNormalVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    peakAcceleration = 0.0;

    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
    expectedDeltaV = 22613.0;
    expectedPeakAcceleration = 56.0e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 100.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );



    /// Fifth Earth-Mercury transfer.

    numberOfRevolutions = 0;

    julianDate = 4355 * physical_constants::JULIAN_DAY;

    timeOfFlight = 160.0;

    // Retrieve cartesian state at departure and arrival.
    cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight  * physical_constants::JULIAN_DAY );

    frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Create base function settings for the components of the radial velocity composite function.
    firstRadialVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdRadialVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear();
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    firstNormalVelocityBaseFunctionSettings = std::make_shared< BaseFunctionHodographicShapingSettings >( );
    secondNormalVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    thirdNormalVelocityBaseFunctionSettings = std::make_shared< PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear();
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( scaledPower, thirdNormalVelocityBaseFunctionSettings ) );


    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    VelocityShapingMethod = shape_based_methods::HodographicShaping(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,    
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 500 );
    peakAcceleration = 0.0;

    for ( int currentStep = 0 ; currentStep <= 500 ; currentStep++ ){
        double currentTime = currentStep * stepSize;

        double currentAcceleration = VelocityShapingMethod.computeCurrentThrustAccelerationVector( currentTime ).norm();
        if ( currentAcceleration > peakAcceleration )
        {
            peakAcceleration = currentAcceleration;
        }

    }


    // Check results consistency w.r.t. thesis from D. Gondelach (ADD PROPER REFERENCE)
    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in Mercury's ephemeris.
    expectedDeltaV = 21766.0;
    expectedPeakAcceleration = 50.3e-4;

    // DeltaV provided with a precision of 1 m/s
    BOOST_CHECK_SMALL( std::fabs(  VelocityShapingMethod.computeDeltaV() - expectedDeltaV ), 100.0 );
    // Peak acceleration provided with a precision 1.0e-6 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakAcceleration - expectedPeakAcceleration ), 1e-5 );

}


//! Test.
BOOST_AUTO_TEST_CASE( test_hodographic_shaping_full_propagation )
{
    using namespace shape_based_methods;

    double numberOfRevolutions = 1.0;

    double julianDate = 2458849.5;

    double timeOfFlight = 500.0;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * physical_constants::JULIAN_DAY );


    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    double scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
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
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsRadialVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsRadialVelocityFunction[ 1 ] = 500.0;

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsNormalVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsNormalVelocityFunction[ 1 ] = - 200.0;

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsAxialVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsAxialVelocityFunction[ 1 ] = 2000.0;

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping VelocityShapingMethod(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY, 1,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );


    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::VectorXd > shapingMethodResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    spice_interface::loadStandardSpiceKernels( );

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );


    // Define body to propagate and central body.
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Vehicle" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Sun" );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    // Define integrator settings.
    double stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 50 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize / 400.0 );

    // Define mass function of the vehicle.
    std::function< double( const double ) > newMassFunction = [ = ]( const double currentTime )
    {
        return 200.0 - 50.0 / ( timeOfFlight * physical_constants::JULIAN_DAY ) * currentTime ;
    };
    bodyMap[ "Vehicle" ]->setBodyMassFunction( newMassFunction );


    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return 200.0;
    };

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >, std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0 );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight * physical_constants::JULIAN_DAY );

    // Create pair of propagator settings (for both forward and backward propagations).
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

    propagatorSettings.first =
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, cartesianStateDepartureBody, terminationConditions.first,
                          propagators::cowell, dependentVariablesToSave );
    propagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, cartesianStateDepartureBody, terminationConditions.second,
              propagators::cowell, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    VelocityShapingMethod.computeShapedTrajectoryAndFullPropagation( bodyMap, specificImpulseFunction, integratorSettings, propagatorSettings,
                                                                      fullPropagationResults, shapingMethodResults, dependentVariablesHistory, false );


    // Check that boundary conditions are still fulfilled when free parameters are added.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin()->second[ i ] - cartesianStateDepartureBody[ i ] )
                           / shapingMethodResults.begin()->second[ i ], 1.0e-8 );
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin()->second[ i ] - cartesianStateArrivalBody[ i ] )
                           / shapingMethodResults.rbegin()->second[ i ], 1.0e-8 );
    }

    // Check results consistency between full propagation and shaped trajectory at departure and arrival.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin()->second[ i ] - fullPropagationResults.begin()->second[ i ] )
                           / shapingMethodResults.begin()->second[ i ], 1.0e-7 );
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin()->second[ i ] - fullPropagationResults.rbegin()->second[ i ] )
                           / shapingMethodResults.rbegin()->second[ i ], 1.0e-7 );
    }

}


//! Test full propagation while propagating the spacecraft mass too.
BOOST_AUTO_TEST_CASE( test_hodographic_shaping_full_propagation_mass_propagation )
{
    using namespace shape_based_methods;

    double numberOfRevolutions = 1.0;

    double julianDate = 2458849.5;

    double timeOfFlight = 500.0;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Retrieve cartesian state at departure and arrival.
    Eigen::Vector6d cartesianStateDepartureBody = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d cartesianStateArrivalBody =
            pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * physical_constants::JULIAN_DAY );


    double frequency = 2.0 * mathematical_constants::PI / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    double scaleFactor = 1.0 / ( timeOfFlight * tudat::physical_constants::JULIAN_DAY );

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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents;
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
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents;
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
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsRadialVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsRadialVelocityFunction[ 1 ] = 500.0;

    // Initialize free coefficients vector for normal velocity function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsNormalVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsNormalVelocityFunction[ 1 ] = - 200.0;

    // Initialize free coefficients vector for axial velocity function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 2 );
    freeCoefficientsAxialVelocityFunction[ 0 ] = 500.0;
    freeCoefficientsAxialVelocityFunction[ 1 ] = 2000.0;

    // Create hodographic-shaping object with defined velocity functions and boundary conditions.
    shape_based_methods::HodographicShaping VelocityShapingMethod(
                cartesianStateDepartureBody,
                cartesianStateArrivalBody,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY, 1,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER,
                radialVelocityFunctionComponents,
                normalVelocityFunctionComponents,
                axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction,
                freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );


    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::VectorXd > shapingMethodResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    spice_interface::loadStandardSpiceKernels( );

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );


    // Define body to propagate and central body.
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Vehicle" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Sun" );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    // Define integrator settings.
    double stepSize = ( timeOfFlight * tudat::physical_constants::JULIAN_DAY ) / static_cast< double >( 50 );
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0, stepSize / 400.0 );

    // Set vehicle mass.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 200.0 );


    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return 3000.0;
    };

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_mass_rate_dependent_variables, "Vehicle" ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >, std::shared_ptr< propagators::PropagationTerminationSettings > >
            terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0 );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight * physical_constants::JULIAN_DAY );

    // Create pair of propagator settings (for both forward and backward propagations).
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > propagatorSettings;

    propagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, cartesianStateDepartureBody, terminationConditions.first,
                          propagators::cowell, dependentVariablesToSave );

    propagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, cartesianStateDepartureBody, terminationConditions.second,
                          propagators::cowell, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    VelocityShapingMethod.computeShapedTrajectoryAndFullPropagation( bodyMap, specificImpulseFunction, integratorSettings, propagatorSettings,
                                                                     fullPropagationResults, shapingMethodResults, dependentVariablesHistory,
                                                                     true );

    // Check that boundary conditions are still fulfilled when free parameters are added.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin()->second[ i ] - cartesianStateDepartureBody[ i ] )
                           / shapingMethodResults.begin()->second[ i ], 1.0e-8 );
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin()->second[ i ] - cartesianStateArrivalBody[ i ] )
                           / shapingMethodResults.rbegin()->second[ i ], 1.0e-8 );
    }

    // Check results consistency between full propagation and shaped trajectory at departure and arrival.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.begin()->second[ i ] - fullPropagationResults.begin()->second[ i ] )
                           / shapingMethodResults.begin()->second[ i ], 1.0e-7 );
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults.rbegin()->second[ i ] - fullPropagationResults.rbegin()->second[ i ] )
                           / shapingMethodResults.rbegin()->second[ i ], 1.0e-7 );
    }

    // Check consistency between current and expected mass rates.
    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin() ; itr != dependentVariablesHistory.end() ; itr++ )
    {
        Eigen::Vector3d currentThrustVector = itr->second.segment( 0, 3 );
        double currentMass = fullPropagationResults.at( itr->first )( 6 );
        double currentMassRate = - itr->second( 3 );
        double expectedMassRate = currentThrustVector.norm() * currentMass /
                ( specificImpulseFunction( itr->first ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

        BOOST_CHECK_SMALL( std::fabs( currentMassRate - expectedMassRate ), 1.0e-15 );

    }


}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
