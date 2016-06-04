#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace unit_tests
{

using namespace orbital_element_conversions;
using namespace unit_conversions;
using namespace reference_frames;

BOOST_AUTO_TEST_SUITE( test_spherical_state_conversions )

//! Test inertial to rotating planetocentric frame transformations, using Matlab script of Erwin
//! Mooij to generate reference data.
BOOST_AUTO_TEST_CASE( testSphericalStateConversions )
{
    // Test case 1: arbitrary rotation
    {
        basic_mathematics::Vector6d cartesianState;
        cartesianState<<-1656517.23153109, -5790058.28764025, -2440584.88186829,
                6526.30784888051, -2661.34558272018, 2377.09572383163;

        double testHeadingAngle = 1.229357188236127;
        double testFlightPathAngle = -0.024894033070522;
        double testLatitude = -0.385027359562548;
        double testLongitude = -1.849449608688977;
        double radius = cartesianState.segment( 0, 3 ).norm( );
        double speed = cartesianState.segment( 3, 3 ).norm( );

        basic_mathematics::Vector6d sphericalOrbitState  = convertCartesianToSphericalOrbitalState(
                    cartesianState );
        basic_mathematics::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitState );
        std::cout<<( sphericalOrbitState.transpose( ) )<<std::endl;

        std::cout<<( reconvertedCartesianState.transpose( ) )<<std::endl<<
                   ( cartesianState ).transpose( )<<std::endl;
        std::cout<<( reconvertedCartesianState.segment( 3, 3 ).norm( ) )<<std::endl<<
                   ( cartesianState ).segment( 3, 3 ).norm( )<<std::endl;
        sleep( 10000.0 );

    }

    // Test case 2: rotation with zero and half pi angles.
    {
        basic_mathematics::Vector6d cartesianState;
        cartesianState<<0.0, 6498098.09700000, 0.0, 0.0, 0.0, 7.438147520000000e+03;

        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = 0.0;
        double testLongitude = mathematical_constants::PI / 2.0;
        double radius = cartesianState.segment( 0, 3 ).norm( );
        double speed = cartesianState.segment( 3, 3 ).norm( );

        basic_mathematics::Vector6d sphericalOrbitState  = convertCartesianToSphericalOrbitalState(
                    cartesianState );
        basic_mathematics::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitState );
        std::cout<<( reconvertedCartesianState - cartesianState ).transpose( )<<std::endl;


    }

    // Test case 3: rotation with zero and half pi angles.
    {
        basic_mathematics::Vector6d cartesianState;
        cartesianState<<0.0, 0.0, 6.498098097000000e3, -7.438147520000000e3, 0.0, 0.0;

        double testHeadingAngle = 0.0;
        double testFlightPathAngle = 0.0;
        double testLatitude = mathematical_constants::PI / 2.0;
        double testLongitude = 0.0;
        double radius = cartesianState.segment( 0, 3 ).norm( );
        double speed = cartesianState.segment( 3, 3 ).norm( );

        basic_mathematics::Vector6d sphericalOrbitState  = convertCartesianToSphericalOrbitalState(
                    cartesianState );
        basic_mathematics::Vector6d reconvertedCartesianState  = convertSphericalOrbitalToCartesianState(
                    sphericalOrbitState );
        std::cout<<( reconvertedCartesianState - cartesianState ).transpose( )<<std::endl;

    }



}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


