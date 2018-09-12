/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelExponentialMapElementConversions.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{
namespace unit_tests
{


//! Test the functionality of the time conversion functions.
BOOST_AUTO_TEST_SUITE( test_USMEM_Element_Conversions )

//! Unit test for conversion Keplerian orbital elements to Unified State Model elements.
BOOST_AUTO_TEST_CASE( testconvertKeplerianToUnifiedStateModelExponentialMapElements )
{
    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    // Setting fraction tolerance for correctness evaluation
    double tolerance = 1.0E-14;

    // Declare gravitational parameter of central body
    const double centralBodyGravitationalParameter = 1.32712440018e20; // [m^3/s^2]

    // Initializing default Keplerian orbit
    Eigen::Vector6d keplerianElements = Eigen::Vector6d::Zero( 6 );
    keplerianElements( semiMajorAxisIndex ) = 1.5e11;
    keplerianElements( eccentricityIndex ) = 0.1;
    keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    keplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

    // Unified state model element vector declaration
    Eigen::Vector7d expectedUnifiedStateModelElements = Eigen::Vector7d::Zero( );
    Eigen::Vector7d computedUnifiedStateModelElements = Eigen::Vector7d::Zero( );

    // Case 1: Elliptical prograde orbit (default case).
    {
        // Default case, so no modification necessary.

        // Expected unified state model elements [m/s,m/s,m/s,-,-,-].
        // (Results obtained by converting quaternions from USM7 unit test to exponential map, with MATLAB code).
        expectedUnifiedStateModelElements( CHodographUSMEMIndex ) = 29894.5892222602;
        expectedUnifiedStateModelElements( Rf1HodographUSMEMIndex ) = -260.548512780222;
        expectedUnifiedStateModelElements( Rf2HodographUSMEMIndex ) = 2978.08312848463;
        expectedUnifiedStateModelElements( e1USMEMIndex ) = -5.13130707826462;
        expectedUnifiedStateModelElements( e2USMEMIndex ) = -0.675549392741421;
        expectedUnifiedStateModelElements( e3USMEMIndex ) = -1.44872034788;
        expectedUnifiedStateModelElements( shadowFlagUSMEMIndex ) = 1.0;

        // Compute unified state model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelExponentialMapElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed unified state model elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );

    }




    // Case 2: Hyperbolic retrograde orbit.
    {
        // Modify Keplerian elements [m,-,rad,rad,rad,rad], i.e. overwrite them.
        keplerianElements( semiMajorAxisIndex ) = -1.5e11;
        keplerianElements( eccentricityIndex ) = 2.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
        keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

        // Set expected unified state model elements [m/s,m/s,m/s,-,-,-].
        // (Results obtained by converting quaternions from USM7 unit test to exponential map, with MATLAB code).
        expectedUnifiedStateModelElements( CHodographUSMEMIndex ) = 17173.1340579794;
        expectedUnifiedStateModelElements( Rf1HodographUSMEMIndex ) = -2993.47450825659;
        expectedUnifiedStateModelElements( Rf2HodographUSMEMIndex ) = 34215.5701963558;
        expectedUnifiedStateModelElements( e1USMEMIndex ) = -3.28605731011794;
        expectedUnifiedStateModelElements( e2USMEMIndex ) = -0.432617652092345;
        expectedUnifiedStateModelElements( e3USMEMIndex ) = -0.0378491401992826;
        expectedUnifiedStateModelElements( shadowFlagUSMEMIndex ) = 1.0;

        // Compute unified state model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelExponentialMapElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed unified state model elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );

    }


    // Case 3: Parabolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( semiLatusRectumIndex ) = 1.5e11;
        keplerianElements( eccentricityIndex ) = 1.0;
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 170.0 );
        keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

        // Set expected unified state model elements [m/s,m/s,m/s,-,-,-].
        // (Results obtained by converting quaternions from USM7 unit test to exponential map, with MATLAB code).
        expectedUnifiedStateModelElements( CHodographUSMEMIndex ) = 29744.7407136119;
        expectedUnifiedStateModelElements( Rf1HodographUSMEMIndex ) = -2592.42496973134;
        expectedUnifiedStateModelElements( Rf2HodographUSMEMIndex ) = 29631.5529950138;
        expectedUnifiedStateModelElements( e1USMEMIndex ) = -0.9433847773697;
        expectedUnifiedStateModelElements( e2USMEMIndex ) = 2.99203425653432;
        expectedUnifiedStateModelElements( e3USMEMIndex ) = -0.27421126568371;
        expectedUnifiedStateModelElements( shadowFlagUSMEMIndex ) = 1.0;

        // Compute unified state model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelExponentialMapElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed unified state model elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );
    }

    // Case 4: Circular prograde orbit with non-zero argument of pericenter, test for error.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.0;
            // Eccentricity is zero, while argument of pericenter is non-zero -> should give error
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );

        // Declare variable indicating whether an exception has been thrown.
        bool isExceptionFound = false;

        // Try computing the unified state model elements and catch the expected runtime error.
        try
        {
            computedUnifiedStateModelElements =
                    convertKeplerianToUnifiedStateModelExponentialMapElements( keplerianElements,
                                                                 centralBodyGravitationalParameter );
        }
        catch( std::runtime_error )
        {
            isExceptionFound = true;
        }

        // Check if runtime error has occured
        BOOST_CHECK( isExceptionFound );
    }

    // Case 5: 0 inclination orbit, test for error because longitude of ascending node is non-zero
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.1;
        keplerianElements( inclinationIndex ) = 0.0;

        // Declare variable indicating whether an exception has been thrown.
        bool isExceptionFound = false;

        // Try computing the unified state model elements and catch the expected runtime error.
        try
        {
            computedUnifiedStateModelElements =
                    convertKeplerianToUnifiedStateModelExponentialMapElements( keplerianElements,
                                                                 centralBodyGravitationalParameter );
        }
        catch( std::runtime_error )
        {
            isExceptionFound = true;
        }

        // Check if runtime error has occured
        BOOST_CHECK( isExceptionFound );
    }

    // Case 6: 180 inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( inclinationIndex ) = PI; // = 180 deg

        // Set expected unified state model elements [m/s,m/s,m/s,-,-,-].
        // (Results obtained by converting quaternions from USM7 unit test to exponential map, with MATLAB code).
        expectedUnifiedStateModelElements( CHodographUSMEMIndex ) = 29894.5892222602;
        expectedUnifiedStateModelElements( Rf1HodographUSMEMIndex ) = -260.548512780222;
        expectedUnifiedStateModelElements( Rf2HodographUSMEMIndex ) = 2978.08312848463;
        expectedUnifiedStateModelElements( e1USMEMIndex ) = -0.944695130614469;
        expectedUnifiedStateModelElements( e2USMEMIndex ) = 2.99619016607469;
        expectedUnifiedStateModelElements( e3USMEMIndex ) = -1.92183978547189e-16;
        expectedUnifiedStateModelElements( shadowFlagUSMEMIndex ) = 1.0;

        // Compute unified state model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelExponentialMapElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Because two elements are near-zero, a close fraction/percentage check will fail.
        // Therefore, 1.0 is added to the elements to avoid this
        expectedUnifiedStateModelElements( e3USMEMIndex ) =
                expectedUnifiedStateModelElements( e3USMEMIndex ) + 1.0;
        computedUnifiedStateModelElements( e3USMEMIndex ) =
                computedUnifiedStateModelElements( e3USMEMIndex ) + 1.0;

        // Check if computed elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );
    }

    // Case 7: 0 eccentricity and inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        keplerianElements( eccentricityIndex ) = 0.0;
        keplerianElements( inclinationIndex ) = 0.0;
        keplerianElements( longitudeOfAscendingNodeIndex ) = 0.0; // Default value because of zero inclination
        keplerianElements( argumentOfPeriapsisIndex ) = 0.0; // Default value because of zero eccentricity

        // Expected unified state model elements [m/s,m/s,m/s,-,-,-].
        // (Results obtained by converting quaternions from USM7 unit test to exponential map, with MATLAB code).
        expectedUnifiedStateModelElements( CHodographUSMEMIndex ) = 29744.7407136119;
        expectedUnifiedStateModelElements( Rf1HodographUSMEMIndex ) = 0;
        expectedUnifiedStateModelElements( Rf2HodographUSMEMIndex ) = 0;
        expectedUnifiedStateModelElements( e1USMEMIndex ) = 0;
        expectedUnifiedStateModelElements( e2USMEMIndex ) = 0;
        expectedUnifiedStateModelElements( e3USMEMIndex ) = 2.96705972839036;
        expectedUnifiedStateModelElements( shadowFlagUSMEMIndex ) = 0.0;

        // Compute unified state model elements.
        computedUnifiedStateModelElements =
                convertKeplerianToUnifiedStateModelExponentialMapElements( keplerianElements,
                                                               centralBodyGravitationalParameter );

        // Check if computed elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedUnifiedStateModelElements,
                                           computedUnifiedStateModelElements, tolerance );
    }

    // Case 8: 200 degree inclination orbit, test for error.
    {
        keplerianElements( inclinationIndex ) = convertDegreesToRadians( 200.0 );
        bool isExceptionFound = false;

        // Try to convert Kepler to unified state model Elements
        try
        {
            computedUnifiedStateModelElements = convertKeplerianToUnifiedStateModelExponentialMapElements
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


//! Unit test for the conversion of unified state model elements to Keplerian elements
BOOST_AUTO_TEST_CASE( testconvertUnifiedStateModelExponentialMapToKeplerianElements )
{
    /* Used procedure:
      Because the Kepler to unified state model elements are verified, a subsequent conversion back
      to Keplerian elements should yield the same outcome as the input Keplerian state. This
      principle is used for verification.
     */

    using namespace orbital_element_conversions;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    // Setting fraction tolerance for correctness evaluation
    double tolerance = 1.0E-14;

    // Declare gravitational parameter of central body
    const double centralBodyGravitationalParameter = 1.32712440018e20; // [m^3/s^2]

    // Initializing default Keplerian orbit
    Eigen::Vector6d expectedKeplerianElements = Eigen::Vector6d::Zero( 6 );
    expectedKeplerianElements( semiMajorAxisIndex ) = 1.5e11;
    expectedKeplerianElements( eccentricityIndex ) = 0.1;
    expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
    expectedKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
    expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    expectedKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 170.0 );

    // Declaring computed output vector.
    Eigen::Vector6d computedKeplerianElements = Eigen::Vector6d::Zero( 6 );



    // Case 1: Elliptical prograde orbit (default case).
    {
        // Default case, so no modification necessary.

        // Convert to unified state model elements and back.
        computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                   centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 2: Hyperbolic retrograde orbit.
    {
        // Modify Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiMajorAxisIndex ) = -1.5e11;
        expectedKeplerianElements( eccentricityIndex ) = 2.0;
        expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 160.0 );
        expectedKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 ); // 170 is above limit

        // Convert to unified state model elements and back.
        computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                   centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 3: Parabolic retrograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiLatusRectumIndex ) = 3.5e11;
        expectedKeplerianElements( eccentricityIndex ) = 1.0;
        expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 90.0 );

        // Convert to unified state model elements and back.
        computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                   centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
       }

    // Case 4: Circular prograde orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiMajorAxisIndex ) = 3.5e11;
        expectedKeplerianElements( eccentricityIndex ) = 0.0;
        expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 70.0 );
        expectedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0; // For e = 0, undefined.

        // Convert to unified state model elements and back.
        computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                   centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 5: 0 inclination orbit,
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( eccentricityIndex ) = 0.3;
        expectedKeplerianElements( inclinationIndex ) = 0.0;
        expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0; // Set to zero as for
        // non-inclined orbit planes, this parameter is undefined

        // Convert to unified state model elements and back.
        computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                   centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 6: 180 inclination orbit, test for error.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiMajorAxisIndex ) = 1.5e15;
        expectedKeplerianElements( inclinationIndex ) = PI;
        expectedKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 240.0 );

        // Declare variable indicating whether an exception has been thrown.
        bool isExceptionFound = false;

        // Try convert to unified state model elements and back and catch the expected runtime error.
        try
        {
            computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                 centralBodyGravitationalParameter ),
                        centralBodyGravitationalParameter );
        }
        catch( std::runtime_error )
        {
            isExceptionFound = true;
        }

        // Check if runtime error has occured
        BOOST_CHECK( isExceptionFound );
    }

    // Case 7: 0 eccentricity and inclination orbit.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( eccentricityIndex ) = 0.0;
            // argument of pericenter was set to 0 in case 4, so no error.
        expectedKeplerianElements( inclinationIndex ) = 0.0;
            // longitude of ascending node was set to 0 in case 5, so no error.

        // Convert to unified state model elements and back
        computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                   centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }

    // Case 8: true anomaly exceeding 180 degrees.
    {
        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        expectedKeplerianElements( semiMajorAxisIndex ) = 1.5e11;
        expectedKeplerianElements( eccentricityIndex ) = 0.1;
        expectedKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 50.0 );
        expectedKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 350.0 );
        expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
        expectedKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 240.0 );

        // Convert to unified state model elements and back
        computedKeplerianElements = convertUnifiedStateModelExponentialMapToKeplerianElements(
                    convertKeplerianToUnifiedStateModelExponentialMapElements( expectedKeplerianElements,
                                                                   centralBodyGravitationalParameter ),
                    centralBodyGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, tolerance );
    }
}
BOOST_AUTO_TEST_SUITE_END( )

} // end namespace unit_tests
} // end namespace tudat
