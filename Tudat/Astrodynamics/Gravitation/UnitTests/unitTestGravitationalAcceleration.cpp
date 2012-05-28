/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
