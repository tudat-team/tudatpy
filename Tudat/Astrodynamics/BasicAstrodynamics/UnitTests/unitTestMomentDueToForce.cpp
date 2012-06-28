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
 *      110622    F.M. Engelen      Creation of code.
 *      110822    D. Dirkx          Removed no longer necessary unit tests.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      120220    D. Dirkx          Moved moment due to force test to separate file after split
 *                                  of moment model base class.
 *      120628    M.I. Ganeff       Boostified unit test.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/momentDueToForceModel.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_MomentDueToForce )

//! Testing moment due to force.
BOOST_AUTO_TEST_CASE( testMomentDueToForceUsingClass )
{
    // Initialize input data.
    const Eigen::Vector3d forceCoefficients( 1.1, 1.2, 1.3 );
    const Eigen::Vector3d momentArm( -2.0, 0.0, 0.0 );
    const double dynamicPressure = 50.0;
    const double referenceArea = 2.2;

    const Eigen::Vector3d intermediateExpectedForce
            = forceCoefficients * dynamicPressure * referenceArea;

    // Compute moment.
    tudat::astrodynamics::moment_models::MomentDueToForceModel momentDueToForce;
    momentDueToForce.computeMoment( intermediateExpectedForce, momentArm );
    const Eigen::Vector3d moment = momentDueToForce.getMomentDueToForce( );

    // Compute expected moment.
    const Eigen::Vector3d expectedMomentDueToForce(
                0.0, -1.0 * momentArm( 0 ) * intermediateExpectedForce( 2 ),
                momentArm( 0 ) * intermediateExpectedForce( 1 ) );

    // Check if computed moment matches expected moment.
    BOOST_CHECK_SMALL( moment( 0 ), std::numeric_limits< double >::min( ) );

    BOOST_CHECK_CLOSE_FRACTION( expectedMomentDueToForce( 1 ), moment( 1 ),
                                std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( expectedMomentDueToForce( 2 ), moment( 2 ),
                                std::numeric_limits< double >::epsilon( ) );
}

//! Testing moment due to force.
BOOST_AUTO_TEST_CASE( testMomentDueToForceUsingFreeFunction )
{
    // Initialize input data.
    const Eigen::Vector3d forceCoefficients( 1.1, 1.2, 1.3 );
    const Eigen::Vector3d momentArm( -2.0, 0.0, 0.0 );
    const double dynamicPressure = 50.0;
    const double referenceArea = 2.2;

    const Eigen::Vector3d intermediateExpectedForce
            = forceCoefficients * dynamicPressure * referenceArea;

    // Compute moment.
    const Eigen::Vector3d moment = tudat::astrodynamics::moment_models::computeMomentDueToForce(
                intermediateExpectedForce, momentArm );

    // Compute expected moment.
    const Eigen::Vector3d expectedMomentDueToForce(
                0.0, -1.0 * momentArm( 0 ) * intermediateExpectedForce( 2 ),
                momentArm( 0 ) * intermediateExpectedForce( 1 ) );

    // Check if computed moment matches expected moment.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedMomentDueToForce, moment,
                                       std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
