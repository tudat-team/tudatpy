/*    Copyright (c) 2010-2016, Delft University of Technology
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
 *      160413    M. Van den Broeck File created
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <tudat/Basics/testMacros.h>

#include <tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelElementConversions.h>

#include <tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>

namespace tudat
{
namespace unit_tests
{



//! Test the functionality of the time conversion functions.
BOOST_AUTO_TEST_SUITE( test_USM_Element_Conversions )

//! Unit test for conversion Keplerian orbital elements to Unified State Model elements.
BOOST_AUTO_TEST_CASE( testConvertKeplerianToUnifiedStateModelElements )
{
    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    // Setting fraction tolerance for correctness evaluation
    double tolerance = 1.0E-14;

    // Declare gravitational parameter of central body
    const double centralBodyGravitationalParameter = 1.32712440018e20; // [m^3/s^2]

    // Initializing default Keplerian orbit
    basic_mathematics::Vector6d keplerianElements = Eigen::VectorXd::Zero( 6 );
    keplerianElements( semiMajorAxisIndex ) = 1.0e7;
    keplerianElements( eccentricityIndex ) = 0.1;
    keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );
    // Unified State Model element vector declaration
    basic_mathematics::Vector6d expectedUnifiedStateModelElements
            = Eigen::VectorXd::Zero( 7 );
    basic_mathematics::Vector6d computedUnifiedStateModelElements
            = Eigen::VectorXd::Zero( 7 );

    // Case 1: Elliptical prograde orbit (default case).
    {
        // Default case, so no modification necessary.

        // Expected Unified State Model elements [m/s,m/s,m/s,-,-,-,-].
        // (Results obtained using code archive B. Rˆmgens (2011)).
        expectedUnifiedStateModelElements( CHodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf1HodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf2HodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon1QuaternionIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon2QuaternionIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon3QuaternionIndex ) =
                basic_mathematics::computeModulo( 6.544984694978736, 2.0 * PI );
        expectedUnifiedStateModelElements( etaQuaternionIndex ) = 0.0;

        // Compute Unified State Model elements.
        basic_mathematics::Vector6d computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );

    }

    // Case 2: Hyperbolic retrograde orbit.
    {
        // Modify Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( semiMajorAxisIndex ) = -1.0e7;
        keplerianElements( eccentricityIndex ) = 2.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
        avoidSingularity = true;
        keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

        // Set expected Unified State Model elements [m/s,m/s,m/s,-,-,-,-]. (Results were calculated by
        // hand).
        expectedUnifiedStateModelElements( CHodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf1HodographIndex )
                = 0.0;
        expectedUnifiedStateModelElements( Rf2HodographIndex )
                = 0.0;
        expectedUnifiedStateModelElements( epsilon1QuaternionIndex )
                = 0.0; //Minor error?
        expectedUnifiedStateModelElements( epsilon2QuaternionIndex )
                = 0.0;
        expectedUnifiedStateModelElements( epsilon3QuaternionIndex )
                = 0.0;
        expectedUnifiedStateModelElements( etaQuaternionIndex )
                = 0.0;

        // Compute Unified State Model elements.
        basic_mathematics::Vector6d computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed Unified State Model elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );

    }

    // Case 3: Parabolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( semiMajorAxisIndex ) = 1.0e7;
        keplerianElements( eccentricityIndex ) = 1.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
        keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

        // Set expected modified equinoctial elements [m/s,m/s,m/s,-,-,-,-]. (Results were calculated by
        // hand)
        expectedUnifiedStateModelElements( CHodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf1HodographIndex )
                = 0.0;
        expectedUnifiedStateModelElements( Rf2HodographIndex )
                = 0.0;
        expectedUnifiedStateModelElements( epsilon1QuaternionIndex )
                = 0.0;
        expectedUnifiedStateModelElements( epsilon2QuaternionIndex )
                = 0.0;
        expectedUnifiedStateModelElements( epsilon3QuaternionIndex )
                = 0.0;
        expectedUnifiedStateModelElements( etaQuaternionIndex ) = 0.0;

        // Compute Unified State Model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed Unified State Model elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );
    }

    // Case 4: Circular prograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );

        // Expected Unified State Model elements [m,-,-,-,-,rad].
        // (Results obtained using code archive B. Rˆmgens (2011)).
        expectedUnifiedStateModelElements( CHodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf1HodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf2HodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon1QuaternionIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon2QuaternionIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon3QuaternionIndex ) =
                basic_mathematics::computeModulo( 9.337511498169663, 2.0 * PI );
        expectedUnifiedStateModelElements( etaQuaternionIndex ) = 0.0;

        // Compute Unified State Model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed Unified State Model elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );

        // Compute Unified State Model elements using direct function.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                             centralBodyGravitationalParameter);

        // Check if computed Unified State Model elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );
    }

    // Case 5: 0 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.1;
        keplerianElements( inclinationIndex ) = 0.0;

        // Expected Unified State Model elements [m/s,m/s,m/s,-,-,-,-].
        // (Results obtained using code archive B. Rˆmgens (2011)).
        expectedUnifiedStateModelElements( CHodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf1HodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf2HodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon1QuaternionIndex ) = 0;
        expectedUnifiedStateModelElements( epsilon2QuaternionIndex ) = 0;
        expectedUnifiedStateModelElements( epsilon3QuaternionIndex ) =
                basic_mathematics::computeModulo( 9.337511498169663, 2.0 * PI );
        expectedUnifiedStateModelElements( etaQuaternionIndex ) = 0.0;

        // Compute Unified State Model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed Unified State Model match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );

    }

    // Case 6: 180 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( inclinationIndex ) = PI; // = 180 deg
        avoidSingularity = true;

        // Set expected Unified State Model elements [m/s,m/s,m/s,-,-,-,-]. (Results were calculated by
        // hand).
        expectedUnifiedStateModelElements( CHodographIndex ) = 0.0;
        expectedUnifiedStateModelElements( Rf1HodographIndex )
                = 0.0;
        expectedUnifiedStateModelElements( Rf2HodographIndex )
                = 0.0;
        expectedUnifiedStateModelElements( epsilon1QuaternionIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon2QuaternionIndex ) = 0.0;
        expectedUnifiedStateModelElements( epsilon3QuaternionIndex )
                = 0.0;
        expectedUnifiedStateModelElements( etaQuaternionIndex ) = 0.0;

        // Compute Unified State Model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        basic_mathematics::Vector6d vectorToAdd
                = ( basic_mathematics::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        basic_mathematics::Vector6d expectedUnifiedStateModelElementsPlusOne =
                expectedUnifiedStateModelElements + vectorToAdd;
        basic_mathematics::Vector6d computedUnifiedStateModelElementsPlusOne =
                computedUnifiedStateModelElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElementsPlusOne,
                                           computedUnifiedStateModelElementsPlusOne, tolerance );

        // Compute Unified State Model elements using direct function.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                             centralBodyGravitationalParameter );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        expectedUnifiedStateModelElementsPlusOne =
                expectedUnifiedStateModelElements + vectorToAdd;
        computedUnifiedStateModelElementsPlusOne =
                computedUnifiedStateModelElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElementsPlusOne,
                                           computedUnifiedStateModelElementsPlusOne, tolerance );
    }

    // Case 7: 0 eccentricity and inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.0;
        keplerianElements( inclinationIndex ) = 0.0;
        avoidSingularity = false;

        // Expected modified equinoctial elements [m,-,-,-,-,rad].
        // (Results obtained using code archive B. Rˆmgens (2011)).
        expectedUnifiedStateModelElements( semiLatusRectumIndex ) = 10000000;
        expectedUnifiedStateModelElements( fElementIndex ) = 0;
        expectedUnifiedStateModelElements( gElementIndex ) = 0;
        expectedUnifiedStateModelElements( hElementIndex ) = 0;
        expectedUnifiedStateModelElements( kElementIndex ) = 0;
        expectedUnifiedStateModelElements( trueLongitudeIndex ) =
                basic_mathematics::computeModulo( 9.337511498169663, 2.0 * PI );

        // Compute modified equinoctial elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        basic_mathematics::Vector6d vectorToAdd
                = ( basic_mathematics::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        basic_mathematics::Vector6d expectedUnifiedStateModelElementsPlusOne =
                expectedUnifiedStateModelElements + vectorToAdd;
        basic_mathematics::Vector6d computedUnifiedStateModelElementsPlusOne =
                computedUnifiedStateModelElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElementsPlusOne,
                                           computedUnifiedStateModelElementsPlusOne, tolerance );

        // Compute modified equinoctial elements using direct function.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelElements( keplerianElements );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        expectedUnifiedStateModelElementsPlusOne =
                expectedUnifiedStateModelElements + vectorToAdd;
        computedUnifiedStateModelElementsPlusOne =
                computedUnifiedStateModelElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElementsPlusOne,
                                           computedUnifiedStateModelElementsPlusOne, tolerance );
    }

    // Case 8: 200 degree inclination orbit, test for error.
    {
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 200.0 );
        bool isExceptionFound = false;

        // Try to calculate retrogradeness
        try
        {
            computedUnifiedStateModelElements = convertKeplerianToUnifiedStateModelElements
                    ( keplerianElements, centralBodyGravitationalParameter );
        }
        // Catch the expected runtime error, and set the boolean flag to true.
        catch ( std::runtime_error )
        {
            isExceptionFound = true;
        }

        // Check value of flag.
        BOOST_CHECK( isExceptionFound );
    }
}

//! Unit test for conversion modified equinoctial elements to Keplerian orbital elements
BOOST_AUTO_TEST_CASE( testConvertUnifiedStateModelToKeplerianElements )
{
    /* Used procedure:
      Because the Kepler to modified equinoctial elements are verified, a subsequent conversion back
      to Keplerian elements should yield the same outcome as the input Keplerian state. This
      principle is used for verification.
     */

    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    // Setting fraction tolerance for correctness evaluation
    double tolerance = 1.0E-14;

    // Initializing default Keplerian orbit
    basic_mathematics::Vector6d expectedKeplerianElements = Eigen::VectorXd::Zero( 6 );
    expectedKeplerianElements( semiMajorAxisIndex ) = 1.0e7;
    expectedKeplerianElements( eccentricityIndex ) = 0.1;
    expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    bool avoidSingularity = false;
    expectedKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    expectedKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

    // Declaring computed output vector.
    basic_mathematics::Vector6d computedKeplerianElements = Eigen::VectorXd::Zero( 6 );

    // Case 1: Elliptical prograde orbit (default case).
    {
        // Default case, so no modification necessary.

        // Convert to modified equinoctial elements and back.
        computedKeplerianElements = convertUnifiedStateModelToKeplerianElements(
                    convertKeplerianToUnifiedStateModelElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 2: Hyperbolic retrograde orbit.
    {
        // Modify Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiMajorAxisIndex ) = -1.0e7;
        expectedKeplerianElements( eccentricityIndex ) = 2.0;
        expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 160.0 );
        avoidSingularity = true;
        expectedKeplerianElements( trueAnomalyIndex )
                = convertDegreesToRadians( 10.0 ); // 170 is above limit

        // Convert to modified equinoctial elements and back.
        computedKeplerianElements = convertUnifiedStateModelToKeplerianElements(
                    convertKeplerianToUnifiedStateModelElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 3: Parabolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiMajorAxisIndex ) = 3.678e7;
        expectedKeplerianElements( eccentricityIndex ) = 1.0;
        expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 90.0 );
        avoidSingularity = true;

        // Convert to modified equinoctial elements and back.
        computedKeplerianElements = convertUnifiedStateModelToKeplerianElements(
                    convertKeplerianToUnifiedStateModelElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 4: Circular prograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( eccentricityIndex ) = 0.0;
        expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 70.0 );
        avoidSingularity = false;
        expectedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0; // For e = 0, undefined.

        // Convert to modified equinoctial elements and back.
        computedKeplerianElements = convertUnifiedStateModelToKeplerianElements(
                    convertKeplerianToUnifiedStateModelElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 5: 0 inclination orbit,
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( eccentricityIndex ) = 0.3;
        expectedKeplerianElements( inclinationIndex ) = 0.0;
        avoidSingularity = false;
        expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0; // Set to zero as for
        // non-inclined orbit planes, this parameter is undefined

        // Convert to modified equinoctial elements and back.
        computedKeplerianElements = convertUnifiedStateModelToKeplerianElements(
                    convertKeplerianToUnifiedStateModelElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 6: 180 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiMajorAxisIndex ) = 1.0e10;
        expectedKeplerianElements( inclinationIndex ) = PI;
        avoidSingularity = true;
        expectedKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 240.0 );

        // Convert to modified equinoctial elements and back.
        computedKeplerianElements = convertUnifiedStateModelToKeplerianElements(
                    convertKeplerianToUnifiedStateModelElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 7: 0 eccentricity and inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( eccentricityIndex ) = 0.0;
        expectedKeplerianElements( inclinationIndex ) = 0.0;
        avoidSingularity = false;

        // Convert to modified equinoctial elements and back
        computedKeplerianElements = convertUnifiedStateModelToKeplerianElements(
                    convertKeplerianToUnifiedStateModelElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }
}

} // namespace unit_tests
} // namespace tudat
