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
 *      110501    L. van der Ham    File created.
 *      110607    L. van der Ham    Make code compatible with Tudat revision 114.
 *      110629    L. van der Ham    Modifications according to comments first code check.
 *      110710    K. Kumar          Restructured code; added subtests.
 *      111024    K. Kumar          Error spotted in L1/L2 tests; locations seem swapped.
 *                                  Tests commented out; needs to be fixed.
 *      111027    K. Kumar          Uncommented out tests as bugs fixed by L. van der Ham.
 *      120307    K. Kumar          Moved file.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120518    K. Kumar          Boostified unit test.
 *      120813    P. Musegaas       Updated unit test to new root finding structure.
 *      130319    K. Kumar          Removed dependence on obsolete central gravity field to new
 *                                  celestial body constants file.
 *
 *    References
 *      Mireles James, J.D. Celestial Mechanics Notes Set 4: The Circular Restricted Three Body
 *          Problem, 2006, http://www.math.utexas.edu/users/jjames/hw4Notes.pdf,
 *          last accessed: 26 May, 2012.
 *      JPL, NASA. Astrodynamic Constants, http://ssd.jpl.nasa.gov/?constants,
 *        last updated: 13 Dec, 2012, last accessed: 19th March, 2013.
 *
 *    Notes
 *      Reference values for position Lagrange libration points are taken from (Mireles James,
 *      2006). There seems to be a bug in the computation of the L3 location! Note that this code
 *      uses the Planet class (in Bodies sub-directory), which is marked for deprecation.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Mathematics/RootFinders/newtonRaphson.h"

namespace tudat
{
namespace unit_tests
{

namespace crtbp = tudat::gravitation::circular_restricted_three_body_problem;
using namespace root_finders;

BOOST_AUTO_TEST_SUITE( test_libration_points )

//! Test if computation of mass parameter is working correctly.
BOOST_AUTO_TEST_CASE( testComputationOfMassParameter )
{
    // Set expected mass parameter for Earth-Moon system.
    const double expectedMassParameter = 0.0121505811805237;

    // Set Earth gravitational parameter.
    const double earthGravitationalParameter
            = tudat::basic_astrodynamics::celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER;

    // Set Moon gravitational parameter (Earth/Moon mass ratio taken from (JPL, 2012).
    const double moonGravitationalParameter = earthGravitationalParameter / 81.30059;

    // Compute mass parameter.
    const double computedMassParameter = crtbp::computeMassParameter(
                earthGravitationalParameter, moonGravitationalParameter );

    // Check if computed value corresponds to expected mass parameter.
    BOOST_CHECK_CLOSE_FRACTION( expectedMassParameter, computedMassParameter, 1.0e-14 );
}

//! Test if computation of location of L1 Lagrange libration point is working correctly.
BOOST_AUTO_TEST_CASE( testComputationOfLocationOfL1LibrationPoint )
{
    // Declare and initialize Earth-Moon mass parameter from (Mireles James, 2006).
    const double earthMoonMassParameter = 0.012277471;

    // Set expected location of L1.
    const Eigen::Vector3d expectedLocationOfL1( 0.83629259089993, 0.0, 0.0 );

    // Declare L1 libration point object with Earth-Moon mass parameter and Newton-Raphson method
    // with 1000 iterations as maximum and 1.0e-14 relative X-tolerance.
    crtbp::LibrationPoint librationPointL1( earthMoonMassParameter,
                                            boost::make_shared< NewtonRaphson >( 1.0e-14, 1000 ) );

    // Compute location of Lagrange libration point.
    librationPointL1.computeLocationOfLibrationPoint( crtbp::LibrationPoint::l1 );

    // Determine location of libration point in Earth-Moon system.
    const Eigen::Vector3d positionOflibrationPointL1
            = librationPointL1.getLocationOfLagrangeLibrationPoint( );

    // Check if computed location of L1 matches expected location.
    BOOST_CHECK_CLOSE_FRACTION( expectedLocationOfL1.x( ),
                                positionOflibrationPointL1.x( ),
                                1.0e-14 );
    BOOST_CHECK_SMALL( positionOflibrationPointL1.y( ), std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( positionOflibrationPointL1.z( ), std::numeric_limits< double >::min( ) );
}

//! Test if computation of location of L2 Lagrange libration point is working correctly.
BOOST_AUTO_TEST_CASE( testComputationOfLocationOfL2LibrationPoint )
{
    // Declare and initialize Earth-Moon mass parameter from (Mireles James, 2006).
    const double earthMoonMassParameter = 0.012277471;

    // Set expected location of L2.
    const Eigen::Vector3d expectedLocationOfL2( 1.15616816590553, 0.0, 0.0 );

    // Declare L2 libration point object with Earth-Moon mass parameter and Newton-Raphson method
    // with 1000 iterations as maximum and 1.0e-14 relative X-tolerance.
    crtbp::LibrationPoint librationPointL2( earthMoonMassParameter,
                                            boost::make_shared< NewtonRaphson >( 1.0e-14, 1000 ) );

    // Compute location of Lagrange libration point.
    librationPointL2.computeLocationOfLibrationPoint( crtbp::LibrationPoint::l2 );

    // Determine location of libration point in Earth-Moon system.
    const Eigen::Vector3d positionOflibrationPointL2
            = librationPointL2.getLocationOfLagrangeLibrationPoint( );

    // Check if computed location of L2 matches expected location.
    BOOST_CHECK_CLOSE_FRACTION( expectedLocationOfL2.x( ),
                                positionOflibrationPointL2.x( ),
                                1.0e-14 );
    BOOST_CHECK_SMALL( positionOflibrationPointL2.y( ), std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( positionOflibrationPointL2.z( ), std::numeric_limits< double >::min( ) );
}

//! Test if computation of location of L3 Lagrange libration point is working correctly.
// THERE IS A BUG IN THIS CASE!
BOOST_AUTO_TEST_CASE( testComputationOfLocationOfL3LibrationPoint )
{
    // Declare and initialize Earth-Moon mass parameter from (Mireles James, 2006).
    const double earthMoonMassParameter = 0.012277471;

    // Set expected location of L3.
    const Eigen::Vector3d expectedLocationOfL3( -1.00511551160689, 0.0, 0.0 );

    // Declare L3 libration point object with Earth-Moon mass parameter and Newton-Raphson method
    // with 1000 iterations as maximum and 1.0e-14 relative X-tolerance.
    crtbp::LibrationPoint librationPointL3( earthMoonMassParameter,
                                            boost::make_shared< NewtonRaphson >( 1.0e-14, 1000 ) );

    // Compute location of Lagrange libration point.
    librationPointL3.computeLocationOfLibrationPoint( crtbp::LibrationPoint::l3 );

    // Determine location of libration point in Earth-Moon system.
    const Eigen::Vector3d positionOflibrationPointL3
            = librationPointL3.getLocationOfLagrangeLibrationPoint( );

    // Check if computed location of L3 matches expected location.
    BOOST_CHECK_CLOSE_FRACTION( expectedLocationOfL3.x( ),
                                positionOflibrationPointL3.x( ),
                                1.0e-2 );
    BOOST_CHECK_SMALL( positionOflibrationPointL3.y( ), std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( positionOflibrationPointL3.z( ), std::numeric_limits< double >::min( ) );
}

//! Test if computation of location of L4 Lagrange libration point is working correctly.
BOOST_AUTO_TEST_CASE( testComputationOfLocationOfL4LibrationPoint )
{
    // Declare and initialize Earth-Moon mass parameter from (Mireles James, 2006).
    const double earthMoonMassParameter = 0.012277471;

    // Set expected location of L4.
    const Eigen::Vector3d expectedLocationOfL4( 0.487722529, 0.86602540378444, 0.0 );

    // Declare L4 libration point object with Earth-Moon mass parameter and Newton-Raphson method
    // with 1000 iterations as maximum and 1.0e-14 relative X-tolerance.
    crtbp::LibrationPoint librationPointL4( earthMoonMassParameter,
                                            boost::make_shared< NewtonRaphson >( 1.0e-14, 1000 ) );

    // Compute location of Lagrange libration point.
    librationPointL4.computeLocationOfLibrationPoint( crtbp::LibrationPoint::l4 );

    // Determine location of libration point in Earth-Moon system.
    const Eigen::Vector3d positionOflibrationPointL4
            = librationPointL4.getLocationOfLagrangeLibrationPoint( );

    // Check if computed location of L4 matches expected location.
    BOOST_CHECK_CLOSE_FRACTION( expectedLocationOfL4.x( ),
                                positionOflibrationPointL4.x( ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( expectedLocationOfL4.y( ),
                                positionOflibrationPointL4.y( ),
                                1.0e-14 );
    BOOST_CHECK_SMALL( positionOflibrationPointL4.z( ), std::numeric_limits< double >::min( ) );
}

//! Test if computation of location of L5 Lagrange libration point is working correctly.
BOOST_AUTO_TEST_CASE( testComputationOfLocationOfL5LibrationPoint )
{
    // Declare and initialize Earth-Moon mass parameter from (Mireles James, 2006).
    const double earthMoonMassParameter = 0.012277471;

    // Set expected location of L5.
    const Eigen::Vector3d expectedLocationOfL5( 0.487722529, -0.86602540378444, 0.0 );

    // Declare L5 libration point object with Earth-Moon mass parameter and Newton-Raphson method
    // with 1000 iterations as maximum and 1.0e-14 relative X-tolerance.
    crtbp::LibrationPoint librationPointL5( earthMoonMassParameter,
                                            boost::make_shared< NewtonRaphson >( 1.0e-14, 1000 ) );

    // Compute location of Lagrange libration point.
    librationPointL5.computeLocationOfLibrationPoint( crtbp::LibrationPoint::l5 );

    // Determine location of libration point in Earth-Moon system.
    const Eigen::Vector3d positionOflibrationPointL5
            = librationPointL5.getLocationOfLagrangeLibrationPoint( );

    // Check if computed location of L5 matches expected location.
    BOOST_CHECK_CLOSE_FRACTION( expectedLocationOfL5.x( ),
                                positionOflibrationPointL5.x( ),
                                1.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( expectedLocationOfL5.y( ),
                                positionOflibrationPointL5.y( ),
                                1.0e-14 );
    BOOST_CHECK_SMALL( positionOflibrationPointL5.z( ), std::numeric_limits< double >::min( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
