/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120226    K. Kumar          File created.
 *
 *    References
 *      Easy calculation. Gravitational Acceleration Tutorial,
 *          http://easycalculation.com/physics/classical-physics
 *          /learn-gravitational-acceleration.php, last accessed: 26th February, 2012.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/gravitationalAccelerationModel.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_gravitational_acceleration )

//! Test if gravitational acceleration is computed correctly.
BOOST_AUTO_TEST_CASE( testGravitationalAcceleration )
{
    // Case 1: Compute gravitational acceleration exerted on surface of Earth
    //         (Easy calculation, 2012).
    {
        // Set gravitational parameter of Earth [m^3 s^-2].
        double gravitationalParameterOfEarth = 6.6726e-11 * 5.9742e24;

        // Set position vector of Earth [m].
        Eigen::Vector3d positionOfEarth = Eigen::Vector3d::Zero( );

        // Set position vector on Earth surface [m].
        Eigen::Vector3d positionOnEarthSurface( 6.3781e6, 0.0, 0.0 );

        // Compute gravitational accelerating acting on Earth's surface [N].
        Eigen::Vector3d gravitationalAccelerationExertedAtEarthSurface
                = astrodynamics::acceleration_models::computeGravitationalAcceleration(
                    positionOnEarthSurface, gravitationalParameterOfEarth, positionOfEarth );

        // Check if computed gravitational force matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( 9.8, gravitationalAccelerationExertedAtEarthSurface.norm( ),
                                    1.0e-4 );
    }

    // Case 2: Compute gravitational acceleration exerted on Lunar surface
    //         (Easy calculation, 2012).
    {
        // Set universal gravitational constant [m^3 kg^-1 s^-2].
        double universalGravitationalConstant = 6.6726e-11;

        // Set mass of Moon [kg].
        double massOfMoon = 7.36e22;

        // Set position vector of Moon [m].
        Eigen::Vector3d positionOfMoon( 12.65, 0.23, -45.78 );

        // Set position vector on surface of Moon [m].
        Eigen::Vector3d positionOfLunarSurface( 0.0, 1735771.89, 0.0 );

        // Compute gravitational accelerating acting on Lunar surface [N].
        Eigen::Vector3d gravitationalAccelerationExertedAtLunarSurface
                = astrodynamics::acceleration_models::computeGravitationalAcceleration(
                    universalGravitationalConstant, positionOfLunarSurface,
                    massOfMoon, positionOfMoon );

        // Check if computed gravitational force matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( 1.63, gravitationalAccelerationExertedAtLunarSurface.norm( ),
                                    1.0e-6 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
