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
 *      B. Rmgens, "Verified Interval Propagation" (2011). MSc thesis,
 *          Delft University of Technology.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{
namespace unit_tests
{

//! Show the functionality of the unit tests.
BOOST_AUTO_TEST_SUITE( test_orbital_element_conversions )

//! Unit test for conversion Keplerian orbital elements to modified equinoctial elements.
BOOST_AUTO_TEST_CASE( testConvertKeplerianToModifiedEquinoctialElements )
{
    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    // Setting fraction tolerance for correctness evaluation
    double tolerance = 1.0E-14;

    // Initializing default Keplerian orbit
    Eigen::Vector6d keplerianElements = Eigen::VectorXd::Zero( 6 );
    keplerianElements( semiMajorAxisIndex ) = 1.0e7;
    keplerianElements( eccentricityIndex ) = 0.1;
    keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    bool avoidSingularity = false;
    keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );
    // Modified equinoctial element vector declaration
    Eigen::Vector6d expectedModifiedEquinoctialElements
            = Eigen::VectorXd::Zero( 6 );
    Eigen::Vector6d computedModifiedEquinoctialElements
            = Eigen::VectorXd::Zero( 6 );

    // Case 1: Elliptical prograde orbit (default case).
    {
        // Default case, so no modification necessary.

        // Expected modified equinoctial elements [m,-,-,-,-,rad].
        // (Results obtained using code archive B. Rmgens (2011)).
        expectedModifiedEquinoctialElements( semiLatusRectumIndex ) = 9900000.0;
        expectedModifiedEquinoctialElements( fElementIndex ) = 0.09961946980917456;
        expectedModifiedEquinoctialElements( gElementIndex ) = 0.008715574274765783;
        expectedModifiedEquinoctialElements( hElementIndex ) = 0.4504186100082874;
        expectedModifiedEquinoctialElements( kElementIndex ) = 0.1206893028076694;
        expectedModifiedEquinoctialElements( trueLongitudeIndex ) =
                basic_mathematics::computeModulo( 6.544984694978736, 2.0 * PI );

        // Compute modified equinoctial elements.
        Eigen::Vector6d computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );

        // Compute modified equinoctial elements using direct function.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );
    }

    // Case 2: Hyperbolic retrograde orbit.
    {
        // Modify Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( semiMajorAxisIndex ) = -1.0e7;
        keplerianElements( eccentricityIndex ) = 2.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
        avoidSingularity = true;
        keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand).
        expectedModifiedEquinoctialElements( semiLatusRectumIndex ) = 3.0e7;
        expectedModifiedEquinoctialElements( fElementIndex )
                = 1.8126155740732999264851053135086;
        expectedModifiedEquinoctialElements( gElementIndex )
                = -0.84523652348139887237395697929546;
        expectedModifiedEquinoctialElements( hElementIndex )
                = 0.0845075596072044152327702959491; //Minor error?
        expectedModifiedEquinoctialElements( kElementIndex )
                = 0.02264373235107538825570191377426;
        expectedModifiedEquinoctialElements( trueLongitudeIndex )
                = 6.0213859193804370403867331512857;

        // Compute modified equinoctial elements.
        Eigen::Vector6d computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );

        // Compute modified equinoctial elements using direct function.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );
    }

    // Case 3: Parabolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( semiMajorAxisIndex ) = 1.0e7;
        keplerianElements( eccentricityIndex ) = 1.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
        avoidSingularity = true;
        keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand)
        expectedModifiedEquinoctialElements( semiLatusRectumIndex ) = 1.0e7;
        expectedModifiedEquinoctialElements( fElementIndex )
                = 0.90630778703664996324255265675432;
        expectedModifiedEquinoctialElements( gElementIndex )
                = -0.42261826174069943618697848964773;
        expectedModifiedEquinoctialElements( hElementIndex )
                = 0.0845075596072044152327702959491;
        expectedModifiedEquinoctialElements( kElementIndex )
                = 0.02264373235107538825570191377426;
        expectedModifiedEquinoctialElements( trueLongitudeIndex )
                = 2.5307274153917778865393516143085;

        // Compute modified equinoctial elements.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );

        // Compute modified equinoctial elements using direct function.
        convertKeplerianToModifiedEquinoctialElements( keplerianElements );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );
    }

    // Case 4: Circular prograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
        avoidSingularity = false;

        // Expected modified equinoctial elements [m,-,-,-,-,rad].
        // (Results obtained using code archive B. Rmgens (2011)).
        expectedModifiedEquinoctialElements( semiLatusRectumIndex ) = 10000000;
        expectedModifiedEquinoctialElements( fElementIndex ) = 0;
        expectedModifiedEquinoctialElements( gElementIndex ) = 0;
        expectedModifiedEquinoctialElements( hElementIndex ) = 0.4504186100082874;
        expectedModifiedEquinoctialElements( kElementIndex ) = 0.1206893028076694;
        expectedModifiedEquinoctialElements( trueLongitudeIndex ) =
                basic_mathematics::computeModulo( 9.337511498169663, 2.0 * PI );

        // Compute modified equinoctial elements.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );

        // Compute modified equinoctial elements using direct function.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );
    }

    // Case 5: 0 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.1;
        keplerianElements( inclinationIndex ) = 0.0;
        avoidSingularity = false;

        // Expected modified equinoctial elements [m,-,-,-,-,rad].
        // (Results obtained using code archive B. Rmgens (2011)).
        expectedModifiedEquinoctialElements( semiLatusRectumIndex ) = 9900000;
        expectedModifiedEquinoctialElements( fElementIndex ) = 0.09961946980917456;
        expectedModifiedEquinoctialElements( gElementIndex ) = 0.008715574274765783;
        expectedModifiedEquinoctialElements( hElementIndex ) = 0;
        expectedModifiedEquinoctialElements( kElementIndex ) = 0;
        expectedModifiedEquinoctialElements( trueLongitudeIndex ) =
                basic_mathematics::computeModulo( 9.337511498169663, 2.0 * PI );

        // Compute modified equinoctial elements.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );

        // Compute modified equinoctial elements using direct function.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements );

        // Check if computed modified equinoctial elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements,
                                           computedModifiedEquinoctialElements, tolerance );
    }

    // Case 6: 180 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( inclinationIndex ) = PI; // = 180 deg
        avoidSingularity = true;

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand).
        expectedModifiedEquinoctialElements( semiLatusRectumIndex ) = 9.9e6;
        expectedModifiedEquinoctialElements( fElementIndex )
                = 0.09063077870366499632425526567543;
        expectedModifiedEquinoctialElements( gElementIndex )
                = -0.04226182617406994361869784896477;
        expectedModifiedEquinoctialElements( hElementIndex ) = 0.0;
        expectedModifiedEquinoctialElements( kElementIndex ) = 0.0;
        expectedModifiedEquinoctialElements( trueLongitudeIndex )
                = 2.5307274153917778865393516143085;

        // Compute modified equinoctial elements.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        Eigen::Vector6d vectorToAdd
                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        Eigen::Vector6d expectedModifiedEquinoctialElementsPlusOne =
                expectedModifiedEquinoctialElements + vectorToAdd;
        Eigen::Vector6d computedModifiedEquinoctialElementsPlusOne =
                computedModifiedEquinoctialElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElementsPlusOne,
                                           computedModifiedEquinoctialElementsPlusOne, tolerance );

        // Compute modified equinoctial elements using direct function.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        expectedModifiedEquinoctialElementsPlusOne =
                expectedModifiedEquinoctialElements + vectorToAdd;
        computedModifiedEquinoctialElementsPlusOne =
                computedModifiedEquinoctialElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElementsPlusOne,
                                           computedModifiedEquinoctialElementsPlusOne, tolerance );
    }

    // Case 7: 0 eccentricity and inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.0;
        keplerianElements( inclinationIndex ) = 0.0;
        avoidSingularity = false;

        // Expected modified equinoctial elements [m,-,-,-,-,rad].
        // (Results obtained using code archive B. Rmgens (2011)).
        expectedModifiedEquinoctialElements( semiLatusRectumIndex ) = 10000000;
        expectedModifiedEquinoctialElements( fElementIndex ) = 0;
        expectedModifiedEquinoctialElements( gElementIndex ) = 0;
        expectedModifiedEquinoctialElements( hElementIndex ) = 0;
        expectedModifiedEquinoctialElements( kElementIndex ) = 0;
        expectedModifiedEquinoctialElements( trueLongitudeIndex ) =
                basic_mathematics::computeModulo( 9.337511498169663, 2.0 * PI );

        // Compute modified equinoctial elements.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                               avoidSingularity );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        Eigen::Vector6d vectorToAdd
                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        Eigen::Vector6d expectedModifiedEquinoctialElementsPlusOne =
                expectedModifiedEquinoctialElements + vectorToAdd;
        Eigen::Vector6d computedModifiedEquinoctialElementsPlusOne =
                computedModifiedEquinoctialElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElementsPlusOne,
                                           computedModifiedEquinoctialElementsPlusOne, tolerance );

        // Compute modified equinoctial elements using direct function.
        computedModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( keplerianElements );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        expectedModifiedEquinoctialElementsPlusOne =
                expectedModifiedEquinoctialElements + vectorToAdd;
        computedModifiedEquinoctialElementsPlusOne =
                computedModifiedEquinoctialElements + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElementsPlusOne,
                                           computedModifiedEquinoctialElementsPlusOne, tolerance );
    }

    // Case 8: 200 degree inclination orbit, test for error.
    {
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 200.0 );
        bool isExceptionFound = false;

        // Try to calculate retrogradeness
        try
        {
            computedModifiedEquinoctialElements = convertKeplerianToModifiedEquinoctialElements
                    ( keplerianElements, avoidSingularity );
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
BOOST_AUTO_TEST_CASE( testConvertModifiedEquinoctialToKeplerianElements )
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
    Eigen::Vector6d expectedKeplerianElements = Eigen::VectorXd::Zero( 6 );
    expectedKeplerianElements( semiMajorAxisIndex ) = 1.0e7;
    expectedKeplerianElements( eccentricityIndex ) = 0.1;
    expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    bool avoidSingularity = false;
    expectedKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    expectedKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

    // Declaring computed output vector.
    Eigen::Vector6d computedKeplerianElements = Eigen::VectorXd::Zero( 6 );

    // Case 1: Elliptical prograde orbit (default case).
    {
        // Default case, so no modification necessary.

        // Convert to modified equinoctial elements and back.
        computedKeplerianElements = convertModifiedEquinoctialToKeplerianElements(
                    convertKeplerianToModifiedEquinoctialElements( expectedKeplerianElements,
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
        computedKeplerianElements = convertModifiedEquinoctialToKeplerianElements(
                    convertKeplerianToModifiedEquinoctialElements( expectedKeplerianElements,
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
        computedKeplerianElements = convertModifiedEquinoctialToKeplerianElements(
                    convertKeplerianToModifiedEquinoctialElements( expectedKeplerianElements,
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
        computedKeplerianElements = convertModifiedEquinoctialToKeplerianElements(
                    convertKeplerianToModifiedEquinoctialElements( expectedKeplerianElements,
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
        computedKeplerianElements = convertModifiedEquinoctialToKeplerianElements(
                    convertKeplerianToModifiedEquinoctialElements( expectedKeplerianElements,
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
        computedKeplerianElements = convertModifiedEquinoctialToKeplerianElements(
                    convertKeplerianToModifiedEquinoctialElements( expectedKeplerianElements,
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
        computedKeplerianElements = convertModifiedEquinoctialToKeplerianElements(
                    convertKeplerianToModifiedEquinoctialElements( expectedKeplerianElements,
                                                                   avoidSingularity ),
                    avoidSingularity );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }
}

//! Unit test for conversion of Cartesian to modified equinoctial elements.
BOOST_AUTO_TEST_CASE( testConvertCartesianElementsToModifiedEquinoctialElements )
{
    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    double tolerance = 1.0E-14;

    Eigen::Vector6d testMEE = Eigen::VectorXd::Zero( 6 );
    Eigen::Vector6d computedMEE = Eigen::VectorXd::Zero( 6 );
    Eigen::Vector6d testCartesianElements = Eigen::VectorXd::Zero( 6 );

    // Set default Keplerian elements [m,-,rad,rad,rad,rad].
    Eigen::Vector6d testKepler = Eigen::VectorXd::Zero( 6 );
    testKepler( semiMajorAxisIndex ) = 1.0e7;
    testKepler( eccentricityIndex ) = 0.1;
    testKepler( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    bool avoidSingularity = false;
    testKepler( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    testKepler( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

    double gravitationalParameter = 398600.44e9; // Earth's, but any parameter would do.

    // Case 1: Elliptical prograde orbit.
    {
        // Default, so no modification necessary

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand)
        testMEE( semiLatusRectumIndex ) = 9.9e6;
        testMEE( fElementIndex )
                = 0.09961946980917455322950104024739;
        testMEE( gElementIndex )
                = 0.00871557427476581735580642708375;
        testMEE( hElementIndex )
                = 0.45041861000828740764931177254188;
        testMEE( kElementIndex )
                = 0.12068930280766941437578622043344;
        testMEE( trueLongitudeIndex )
                = 3.0543261909900767596164588448551;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                     gravitationalParameter );

        // Convert to modified equinoctial elements.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter,
                                                                     avoidSingularity );

        // Compare, because element 2 is quite small, tolerance is less stringent than usual.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, 1.0E-13 );

        // Convert to modified equinoctial elements using direct function.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter );

        // Compare, because element 2 is quite small, tolerance is less stringent than usual.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, 1.0E-13 );
    }

    // Case 2: Hyperbolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( semiMajorAxisIndex ) = -1.0e7;
        testKepler( eccentricityIndex ) = 2.0;
        testKepler( inclinationIndex ) = convertDegreesToRadians( 170.0 );
        avoidSingularity = true;
        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand)
        testMEE( semiLatusRectumIndex ) = 3.0e7;
        testMEE( fElementIndex )
                = 1.8126155740732999264851053135086;
        testMEE( gElementIndex )
                = -0.84523652348139887237395697929546;
        testMEE( hElementIndex )
                = 0.0845075596072044152327702959491;
        testMEE( kElementIndex )
                = 0.02264373235107538825570191377426;
        testMEE( trueLongitudeIndex )
                = 6.0213859193804370403867331512857;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                     gravitationalParameter );

        // Convert to modified equinoctial elements.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter,
                                                                     avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );

        // Convert to modified equinoctial elements using direct function.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );
    }

    // Case 3: Parabolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( semiMajorAxisIndex ) = 1.0e7;
        testKepler( eccentricityIndex ) = 1.0;
        avoidSingularity = true;
        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand).
        testMEE( semiLatusRectumIndex ) = 1.0e7;
        testMEE( fElementIndex )
                = 0.90630778703664996324255265675432;
        testMEE( gElementIndex )
                = -0.42261826174069943618697848964773;
        testMEE( hElementIndex )
                = 0.0845075596072044152327702959491;
        testMEE( kElementIndex )
                = 0.02264373235107538825570191377426;
        testMEE( trueLongitudeIndex )
                = 2.5307274153917778865393516143085;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                     gravitationalParameter );

        // Convert to modified equinoctial elements.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter,
                                                                     avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );

        // Convert to modified equinoctial elements using direct function.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testMEE, computedMEE, tolerance );
    }

    // Case 4: Circular prograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler ( eccentricityIndex ) = 0.0;
        testKepler ( inclinationIndex ) = convertDegreesToRadians( 50.0 );
        avoidSingularity = false;
        testKepler ( argumentOfPeriapsisIndex ) = 0.0; // e = 0, so actually undefined

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand).
        testMEE ( semiLatusRectumIndex ) = 1.0e7;
        testMEE ( fElementIndex ) = 0.0;
        testMEE ( gElementIndex ) = 0.0;
        testMEE ( hElementIndex )
                = 0.45041861000828740764931177254188;
        testMEE ( kElementIndex )
                = 0.12068930280766941437578622043344;
        testMEE ( trueLongitudeIndex )
                = 3.2288591161895097173088279217039;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                     gravitationalParameter );

        // Convert to modified equinoctial elements.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter,
                                                                     avoidSingularity );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        Eigen::Vector6d vectorToAdd
                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );

        // Convert to modified equinoctial elements using direct function
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this.
        computedMeePlusOne = computedMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
    }

    // Case 5: 0 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler ( eccentricityIndex ) = 0.1;
        testKepler ( inclinationIndex ) = 0.0;
        avoidSingularity = false;
        testKepler ( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 260.0 );
        testKepler ( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand).
        testMEE( semiLatusRectumIndex ) = 9.9e6;
        testMEE( fElementIndex )
                = -0.01736481776669303488517166267693;
        testMEE( gElementIndex )
                = -0.09848077530122080593667430245895;
        testMEE( hElementIndex ) = 0.0;
        testMEE( kElementIndex ) = 0.0;
        testMEE( trueLongitudeIndex )
                = 1.221730476396030703846583537942;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                     gravitationalParameter );

        // Convert to modified equinoctial elements.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter,
                                                                     avoidSingularity );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this.
        Eigen::Vector6d vectorToAdd
                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );

        // Convert to modified equinoctial elements using direct function.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this.
        computedMeePlusOne = computedMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
    }

    // Case 6: 180 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( eccentricityIndex ) = 0.1;
        testKepler( inclinationIndex ) = PI; // = 180 deg
        avoidSingularity = true;
        testKepler( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 12.0 );
        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 190.0 );

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand).
        testMEE ( semiLatusRectumIndex ) = 9.9e6;
        testMEE ( fElementIndex )
                = 0.09781476007338056379285667478696;
        testMEE ( gElementIndex )
                = 0.02079116908177593371017422844051;
        testMEE ( hElementIndex )
                = 0.0;
        testMEE ( kElementIndex )
                = 0.0;
        testMEE ( trueLongitudeIndex )
                = 3.525565089028545745385855352347;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                     gravitationalParameter );

        // Convert to modified equinoctial elements.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter,
                                                                     avoidSingularity );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this.
        Eigen::Vector6d vectorToAdd
                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );

        // Convert to modified equinoctial elements using direct function.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this.
        computedMeePlusOne = computedMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
    }

    // Case 7: 0 eccentricity and inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( eccentricityIndex ) = 0.0;
        testKepler( inclinationIndex ) = 0.0;
        avoidSingularity = false;

        // Set expected modified equinoctial elements [m,-,-,-,-,rad]. (Results were calculated by
        // hand).
        testMEE ( semiLatusRectumIndex ) = 1.0e7; // Circular
        testMEE ( fElementIndex ) = 0.0;
        testMEE ( gElementIndex ) = 0.0;
        testMEE ( hElementIndex ) = 0.0;
        testMEE ( kElementIndex ) = 0.0;
        testMEE ( trueLongitudeIndex )
                = 3.525565089028545745385855352347;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        testCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                     gravitationalParameter );

        // Convert to modified equinoctial elements.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter,
                                                                     avoidSingularity );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this.
        Eigen::Vector6d vectorToAdd
                = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished( );
        Eigen::Vector6d computedMeePlusOne = computedMEE + vectorToAdd;
        Eigen::Vector6d testMeePlusOne = testMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );

        // Convert to modified equinoctial elements using direct function.
        computedMEE = convertCartesianToModifiedEquinoctialElements( testCartesianElements,
                                                                     gravitationalParameter );

        // Check if computed elements match the expected values.
        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this.
        computedMeePlusOne = computedMEE + vectorToAdd;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedMeePlusOne, testMeePlusOne, tolerance );
    }
}

//! Unit test for conversion of modified equinoctial elements to Cartesian.
BOOST_AUTO_TEST_CASE( testConvertModifiedEquinoctialToCartesianElements )
{
    /* Used procedure:
      The Cartesian expected outcome is computed from the verified Kepler to Cartesian conversion.
      Subsequently, the Kepler state is converted to modified equinoctial elements and then
      converted back to Cartesian elements. Outcomes are compared.
     */

    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    // Tolerance precision: two orders higher than machine due to three conversions being applied
    // (accumulation of error) in order to save on manual labor.
    double tolerance = 2.0E-14;

    Eigen::Vector6d intermediateModifiedEquinoctialElements
            = Eigen::VectorXd::Zero( 6 );
    Eigen::Vector6d expectedCartesianElements = Eigen::VectorXd::Zero( 6 );
    Eigen::Vector6d computedCartesianElements = Eigen::VectorXd::Zero( 6 );

    // Set default Keplerian elements [m,-,rad,rad,rad,rad].
    Eigen::Vector6d testKepler = Eigen::VectorXd::Zero( 6 );
    testKepler( semiMajorAxisIndex ) = 1.0e7;
    testKepler( eccentricityIndex ) = 0.1;
    testKepler( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    bool avoidSingularity = false;
    testKepler( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    testKepler( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

    double gravitationalParameter = 398600.44e9; // Earth's, but any parameter would do.

    // Case 1: Elliptical prograde orbit.
    {
        // Default, so no modification necessary.

        // Create expected Cartesian vector through the verified Kepler to Cartesian routine.
        expectedCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                         gravitationalParameter );

        // Convert to modified equinoctial elements, then that to Cartesian.
        intermediateModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( testKepler, avoidSingularity );
        computedCartesianElements = convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );

        computedCartesianElements = convertModifiedEquinoctialToCartesianElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );
    }

    // Case 2: Hyperbolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( semiMajorAxisIndex ) = -1.0e7;
        testKepler( eccentricityIndex ) = 2.0;
        testKepler( inclinationIndex )
                = convertDegreesToRadians( 170.0 ); // Between 90 and 180 is retrograde
        avoidSingularity = true;
        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        expectedCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                         gravitationalParameter );

        // Convert to modified equinoctial elements, then to Cartesian.
        intermediateModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( testKepler, avoidSingularity );
        computedCartesianElements = convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Check if computed elements match the expected values.
        // Because element z is ~10^16 smaller than the other elements, it is only checked whether
        // its value is 'sufficiently' close to zero.
        BOOST_CHECK_SMALL( expectedCartesianElements( 2 ), 1.0E-9 );
        BOOST_CHECK_SMALL( computedCartesianElements( 2 ), 1.0E-9 );
        expectedCartesianElements( 2 ) = 0.0;
        computedCartesianElements( 2 ) = 0.0;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );
    }

    // Case 3: Parabolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( semiMajorAxisIndex ) = 1.0e7;
        testKepler( eccentricityIndex ) = 1.0;
        avoidSingularity = true;
        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        expectedCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                         gravitationalParameter );

        // Convert to modified equinoctial elements, then to Cartesian.
        intermediateModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( testKepler, avoidSingularity );
        computedCartesianElements = convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );
    }

    // Case 4: Circular prograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler ( eccentricityIndex ) = 0.0;
        testKepler ( inclinationIndex ) = convertDegreesToRadians( 50.0 );
        avoidSingularity = false;
        testKepler ( argumentOfPeriapsisIndex ) = 0.0; // e = 0, so actually undefined.

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        expectedCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                         gravitationalParameter );

        // Convert to modified equinoctial elements, then to Cartesian.
        intermediateModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( testKepler, avoidSingularity );
        computedCartesianElements = convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           1.0E-13 );

        computedCartesianElements = convertModifiedEquinoctialToCartesianElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           1.0E-13 );
    }

    // Case 5: 0 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler ( eccentricityIndex ) = 0.1;
        testKepler ( inclinationIndex ) = 0.0;
        avoidSingularity = false;
        testKepler ( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 260.0 );
        testKepler ( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 0.0 );

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        expectedCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                         gravitationalParameter );

        // Convert to modified equinoctial elements, then to Cartesian.
        intermediateModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( testKepler, avoidSingularity );
        computedCartesianElements = convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );

        computedCartesianElements = convertModifiedEquinoctialToCartesianElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );
    }

    // Case 6: 180 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( eccentricityIndex ) = 0.1;
        testKepler( inclinationIndex ) = PI; // = 180 deg
        avoidSingularity = true;
        testKepler( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 12.0 );
        testKepler( trueAnomalyIndex ) = convertDegreesToRadians( 190.0 );

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        expectedCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                         gravitationalParameter );

        // Convert to modified equinoctial elements, then to Cartesian.
        intermediateModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( testKepler, avoidSingularity );
        computedCartesianElements = convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );

//        computedCartesianElements = convertModifiedEquinoctialToCartesianElements(
//                    intermediateModifiedEquinoctialElements, gravitationalParameter,
//                    avoidSingularity );

//        // Compare.
//        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
//                                           tolerance );
    }

    // Case 7: 0 eccentricity and inclination.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        testKepler( eccentricityIndex ) = 0.0;
        testKepler( inclinationIndex ) = 0.0;
        avoidSingularity = false;

        // Create starting Cartesian vector through the verified Kepler to Cartesian routine.
        expectedCartesianElements = convertKeplerianToCartesianElements( testKepler,
                                                                         gravitationalParameter );

        // Convert to modified equinoctial elements, then to Cartesian.
        intermediateModifiedEquinoctialElements =
                convertKeplerianToModifiedEquinoctialElements( testKepler, avoidSingularity );
        computedCartesianElements = convertModifiedEquinoctialToCartesianElementsViaKeplerElements(
                    intermediateModifiedEquinoctialElements, gravitationalParameter,
                    avoidSingularity );

        // Compare.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements, computedCartesianElements,
                                           tolerance );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
