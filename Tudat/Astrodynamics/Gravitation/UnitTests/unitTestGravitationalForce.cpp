/*    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      120209    K. Kumar          File created.
 *
 *    References
 *      Easy calculation. Newton's Law of Gravity Tutorial,
 *          http://easycalculation.com/physics/classical-physics/learn-newtons-law.php, last
 *          accessed: 12th February, 2012.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/gravitationalForceModel.h"

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
                = tudat::force_models::computeGravitationalForce(
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
                = tudat::force_models::computeGravitationalForce(
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
