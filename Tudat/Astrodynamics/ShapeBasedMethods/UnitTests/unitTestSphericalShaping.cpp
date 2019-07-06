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

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/sphericalShaping.h"

namespace tudat
{
namespace unit_tests
{

//! Test spherical shaping implementation.
BOOST_AUTO_TEST_SUITE( test_spherical_shaping )

//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping1 )
{

    double numberOfRevolutions = 1.0;
    double julianDate = 8174.5 * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    std::cout << "cartesian state departure body = Earth: " << pointerToDepartureBodyEphemeris->getCartesianState( julianDate  ) << "\n\n";
    std::cout << "cartesian state arrival body = Mars: "
            << pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * physical_constants::JULIAN_DAY ) << "\n\n";

    Eigen::VectorXd radialFunctionCoefficients = ( Eigen::Vector7d() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished();
    Eigen::Vector4d elevationFunctionCoefficients = ( Eigen::Vector4d() << 1.0, 1.0, 1.0, 1.0 ).finished();

    Eigen::Vector6d normalisedInitialState;
    normalisedInitialState.segment( 0, 3 ) = pointerToDepartureBodyEphemeris->getCartesianState( julianDate ).segment( 0, 3 )
            / physical_constants::ASTRONOMICAL_UNIT;
    normalisedInitialState.segment( 3, 3 ) = pointerToDepartureBodyEphemeris->getCartesianState( julianDate ).segment( 3, 3 )
             * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    Eigen::Vector6d normalisedFinalState;
    normalisedFinalState.segment( 0, 3 ) = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ).segment( 0, 3 )
            / physical_constants::ASTRONOMICAL_UNIT;
    normalisedFinalState.segment( 3, 3 ) = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ).segment( 3, 3 )
             * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                normalisedInitialState,
                normalisedFinalState,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY  / physical_constants::JULIAN_YEAR, 0.000703 /*- 0.00000000000007*/,
                radialFunctionCoefficients, elevationFunctionCoefficients,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER  * std::pow( physical_constants::JULIAN_YEAR, 2.0 )
                / std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ),
                root_finders::bisection_root_finder, 1.0e-6, 1.0e-1, TUDAT_NAN, 30, 1.0e-6 );

    std::cout << "test initial spherical position: " << sphericalShaping.computeCurrentSphericalState( sphericalShaping.getInitialAzimuthalAngle() /*2.25776*/ ).transpose() << "\n\n";
    std::cout << "test final spherical position: " << sphericalShaping.computeCurrentSphericalState( sphericalShaping.getFinalAzimuthalAngle() /*12.8632*/ ).transpose() << "\n\n";

    std::cout << "test initial cartesian position: " << sphericalShaping.computeCurrentCartesianState( sphericalShaping.getInitialAzimuthalAngle() /*2.25776*/ ).transpose() << "\n\n";
    std::cout << "test final cartesian position: " << sphericalShaping.computeCurrentCartesianState( sphericalShaping.getFinalAzimuthalAngle() /*12.8632*/ ).transpose() << "\n\n";

    std::cout << "difference final cartesian state: " << ( sphericalShaping.computeCurrentCartesianState( sphericalShaping.getFinalAzimuthalAngle() )
                 - pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ) ).transpose() << "\n\n";

    std::map< double, Eigen::Vector6d > mapPosition;
    std::map< double, Eigen::Vector3d > mapAcceleration;

    // Compute step size.
    double stepSize = ( sphericalShaping.getFinalAzimuthalAngle() - sphericalShaping.getInitialAzimuthalAngle() ) / 5000.0;

    // Initialise peak acceleration.
    double peakThrustAcceleration = 0.0;

    // Check that the trajectory is feasible, ie curved toward the central body.
    for ( int i = 0 ; i <= 5000 ; i++ )
    {
        double currentThetaAngle = sphericalShaping.getInitialAzimuthalAngle() + i * stepSize;

        Eigen::Vector6d currentUnnormalisedState;
        currentUnnormalisedState.segment( 0, 3 ) = sphericalShaping.computeCurrentCartesianState( currentThetaAngle ).segment( 0, 3 )
                * physical_constants::ASTRONOMICAL_UNIT;
        currentUnnormalisedState.segment( 3, 3 ) = sphericalShaping.computeCurrentCartesianState( currentThetaAngle ).segment( 3, 3 )
                * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;

        mapPosition[ currentThetaAngle ] = currentUnnormalisedState;
        mapAcceleration[ currentThetaAngle ] = sphericalShaping.computeSphericalControlAccelerationVector( currentThetaAngle )
                * physical_constants::ASTRONOMICAL_UNIT / std::pow( physical_constants::JULIAN_YEAR, 2.0 );

        if ( sphericalShaping.computeSphericalControlAccelerationVector( currentThetaAngle ).norm() >  peakThrustAcceleration )
        {
            peakThrustAcceleration = sphericalShaping.computeSphericalControlAccelerationVector( currentThetaAngle ).norm();
        }
    }

    tudat::input_output::writeDataMapToTextFile( mapPosition,
                                          "mapPosition.dat",
                                           "C:/Users/chamb/Documents/Master_2/SOCIS/SphericalShapingTest/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    tudat::input_output::writeDataMapToTextFile( mapAcceleration,
                                          "mapAcceleration.dat",
                                           "C:/Users/chamb/Documents/Master_2/SOCIS/SphericalShapingTest/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    std::cout << "deltaV: " << sphericalShaping.computeDeltav() * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR << "\n\n";
    std::cout << "peak acceleration: " << peakThrustAcceleration * physical_constants::ASTRONOMICAL_UNIT / std::pow( physical_constants::JULIAN_YEAR, 2.0 ) << "\n\n";
    std::cout << "free coefficient: " << sphericalShaping.getRadialCompositionFunctionCoefficients().transpose() << "\n\n";

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
