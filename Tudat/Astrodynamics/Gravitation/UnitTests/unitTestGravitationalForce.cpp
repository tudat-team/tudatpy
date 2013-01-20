/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      120209    K. Kumar          File created.
 *
 *    References
 *      Easy calculation. Newton's Law of Gravity Tutorial,
 *          http://easycalculation.com/physics/classical-physics/learn-newtons-law.php, last
 *          accessed: 12th February, 2012.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_gravitational_force )

//! Test if gravitational force is computed correctly.
BOOST_AUTO_TEST_CASE( testGravitationalForce )
{
    // Case 1: Compute gravitational force exerted on boy, due to Earth (Easy calculation, 2012).
    {
        // Set gravitational parameter of Earth [m^3 s^-2].
        double gravitationalParameterOfEarth = 6.6726e-11 * 5.98e24;

        // Set position vector of Earth [m].
        Eigen::Vector3d positionOfEarth = Eigen::Vector3d::Zero( );

        // Set mass of boy [kg]
        double massOfBoy = 70.0;

        // Set position vector of boy [m].
        Eigen::Vector3d positionOfBoy( 6.38e6, 0.0, 0.0 );

        // Compute gravitational force acting on boy [N].
        Eigen::Vector3d gravitationalForceExertedOnBoy
                = gravitation::computeGravitationalForce(
                    massOfBoy, positionOfBoy,
                    gravitationalParameterOfEarth, positionOfEarth );

        // Check if computed gravitational force matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( 685.54, gravitationalForceExertedOnBoy.norm( ), 1.0e-3 );
    }

    // Case 2: Compute gravitational force exerted on arbitrary body1, due to arbitrary body2
    //         (Easy calculation, 2012).
    {
        // Set universal gravitational constant [m^3 kg^-1 s^-2].
        double universalGravitationalConstant = 6.6726e-11;

        // Set mass of body1 [kg].
        double massOfBody1 = 1.0e4;

        // Set position vector of body1 [m].
        Eigen::Vector3d positionOfBody1( 0.0, 5.0, 0.0 );

        // Set mass of body2 [kg].
        double massOfBody2 = 2.0e4;

        // Set position vector of body2 [m].
        Eigen::Vector3d positionOfBody2( 0.0, -5.0, 0.0 );

        // Compute gravitational force acting on body1 [N].
        Eigen::Vector3d gravitationalForceExertedOnBody1
                = gravitation::computeGravitationalForce(
                    universalGravitationalConstant, massOfBody1,
                    positionOfBody1, massOfBody2, positionOfBody2 );

        // Check if computed gravitational force matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( 13345200.0e-11, gravitationalForceExertedOnBody1.norm( ),
                                    std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
