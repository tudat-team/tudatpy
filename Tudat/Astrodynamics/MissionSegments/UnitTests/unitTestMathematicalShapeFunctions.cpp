/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The current test cases are used to check the code of the mathematical functions, each case
 *      represents a special case where certain elements of the function can be eliminated.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/MissionSegments/improvedInversePolynomialWall.h"
#include "Tudat/Astrodynamics/MissionSegments/oscillatingFunctionNovak.h"

namespace tudat
{
namespace unit_tests
{

//! Test the spherical shape function code.
BOOST_AUTO_TEST_SUITE( test_spherical_shape_function )

//! Test case: inverse polynomial function.
BOOST_AUTO_TEST_CASE( ImprovedInversePolynomialWall )
{
    // Using declaration.
    using mathematical_constants::PI;

    // Initialize parameters of the inverse polynomial function.
    std::pair< Eigen::Vector3d, Eigen::Vector3d > inversePolynomialParameters;
    double timeDepParameter;
    double azimuthalAngle = 0.0;

    // Error tolerance.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    //*********************************************************************************************
    // CASE 1: b = 1.0, polar angle = 0.0, all other parameters 0.0
    //*********************************************************************************************

    // Case specific changes.
    inversePolynomialParameters.first = Eigen::Vector3d::Zero(  );
    inversePolynomialParameters.second = Eigen::Vector3d::Zero(  );
    timeDepParameter = 0.0;
    inversePolynomialParameters.first( 1 ) = 1.0; // b

    // Expected results.
    double expectedFunctionValue = 1.0;
    double expectedFirstDerivative = 0.0;
    double expectedSecondDerivative = 1.0;
    double expectedThirdDerivative = 0.0;

    // Initialize mathematical function.
    mission_segments::ImprovedInversePolynomialWall myFunctionCase1(
                boost::lambda::constant( timeDepParameter ) ,
                boost::lambda::constant( inversePolynomialParameters ) );

    // Calculate function value and derivatives.
    double functionValue = myFunctionCase1.evaluate( azimuthalAngle );
    double firstDerivative = myFunctionCase1.computeDerivative( 1 , azimuthalAngle );
    double secondDerivative = myFunctionCase1.computeDerivative( 2 , azimuthalAngle );
    double thirdDerivative = myFunctionCase1.computeDerivative( 3 , azimuthalAngle );

    // Test if the computed values correspond to the expected values, within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( functionValue,
                                expectedFunctionValue,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( firstDerivative,
                                expectedFirstDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( secondDerivative,
                                expectedSecondDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( thirdDerivative,
                                expectedThirdDerivative,
                                tolerance );

    //*********************************************************************************************
    // CASE 2: c = 0.0, polar angle = 0.0, all other parameters 1.0
    //*********************************************************************************************

    // Case-specific changes.
    inversePolynomialParameters.first = Eigen::Vector3d::Ones(  );
    inversePolynomialParameters.second = Eigen::Vector3d::Ones(  );
    timeDepParameter = 1.0;
    inversePolynomialParameters.first( 2 ) = 0.0; // c

    // Expected results.
    expectedFunctionValue = 1.0 / 2.0;
    expectedFirstDerivative = 0.0;
    expectedSecondDerivative = 1.0 / 4.0;
    expectedThirdDerivative = - 6.0 / 4.0;

    // Initialize mathematical function.
    mission_segments::ImprovedInversePolynomialWall myFunctionCase2(
                boost::lambda::constant( timeDepParameter ) ,
                boost::lambda::constant( inversePolynomialParameters ) );

    // Calculate function value and derivatives.
    functionValue = myFunctionCase2.evaluate( azimuthalAngle );
    firstDerivative = myFunctionCase2.computeDerivative( 1 , azimuthalAngle );
    secondDerivative = myFunctionCase2.computeDerivative( 2 , azimuthalAngle );
    thirdDerivative = myFunctionCase2.computeDerivative( 3 , azimuthalAngle );

    // Test if the computed values correspond to the expected values, within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( functionValue,
                                expectedFunctionValue,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( firstDerivative,
                                expectedFirstDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( secondDerivative,
                                expectedSecondDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( thirdDerivative,
                                expectedThirdDerivative,
                                tolerance );

    //*********************************************************************************************
    // CASE 3: c = pi - 1.0, polar angle = 1.0, all other parameters 1.0
    //*********************************************************************************************

    // Case-specific changes.
    inversePolynomialParameters.first( 2 ) = PI - 1.0; // c
    azimuthalAngle = 1.0;

    // Expected results.
    expectedFunctionValue = 1.0 / 4.0;
    expectedFirstDerivative = -9.0 / 8.0;
    expectedSecondDerivative = 93.0 / 16.0;
    expectedThirdDerivative = -267.0 / 8.0;

    // Initialize mathematical function.
    mission_segments::ImprovedInversePolynomialWall myFunctionCase3(
                boost::lambda::constant( timeDepParameter ) ,
                boost::lambda::constant( inversePolynomialParameters ) );

    // Calculate function value and derivatives.
    functionValue = myFunctionCase3.evaluate( azimuthalAngle );
    firstDerivative = myFunctionCase3.computeDerivative( 1 , azimuthalAngle );
    secondDerivative = myFunctionCase3.computeDerivative( 2 , azimuthalAngle );
    thirdDerivative = myFunctionCase3.computeDerivative( 3 , azimuthalAngle );

    // Test if the computed values correspond to the expected values, within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( functionValue,
                                expectedFunctionValue,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( firstDerivative,
                                expectedFirstDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( secondDerivative,
                                expectedSecondDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( thirdDerivative,
                                expectedThirdDerivative,
                                tolerance );
}

//! Test case: oscillating function.
BOOST_AUTO_TEST_CASE( OscillatingFunctionNovak )
{
    // Using declaration.
    using mathematical_constants::PI;

    // Initialize parameters of the inverse polynomial function.
    std::pair< Eigen::Vector2d, Eigen::Vector2d > oscillatingFunctionParameters;
    double azimuthalAngle = 0.0;

    // Error tolerance.
    const double tolerance = std::numeric_limits< double >::epsilon(  );

    //*********************************************************************************************
    // CASE 1: theta = 0.0, all parameters = 1.0
    //*********************************************************************************************

    // Case specific changes.
    oscillatingFunctionParameters.first = Eigen::Vector2d::Ones(  );
    oscillatingFunctionParameters.second = Eigen::Vector2d::Ones(  );

    // Expected results.
    double expectedFunctionValue = 1.0;
    double expectedFirstDerivative = 2.0;
    double expectedSecondDerivative = 1.0;
    double expectedThirdDerivative = - 4.0;

    // Initialize mathematical function.
    mission_segments::OscillatingFunctionNovak myFunctionCase1(
                boost::lambda::constant( oscillatingFunctionParameters ) );

    // Calculate function value and derivatives.
    double functionValue = myFunctionCase1.evaluate( azimuthalAngle );
    double firstDerivative = myFunctionCase1.computeDerivative( 1 , azimuthalAngle );
    double secondDerivative = myFunctionCase1.computeDerivative( 2 , azimuthalAngle );
    double thirdDerivative = myFunctionCase1.computeDerivative( 3 , azimuthalAngle );

    // Test if the computed values correspond to the expected values, within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( functionValue,
                                expectedFunctionValue,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( firstDerivative,
                                expectedFirstDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( secondDerivative,
                                expectedSecondDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( thirdDerivative,
                                expectedThirdDerivative,
                                tolerance );

    //*********************************************************************************************
    // CASE 2: theta = pi/2, all parameters = 1.0
    //*********************************************************************************************

    // Case-specific changes.
    azimuthalAngle = PI / 2.0;

    // Expected results.
    expectedFunctionValue = 1.0 + PI / 2.0;
    expectedFirstDerivative = -PI / 2.0;
    expectedSecondDerivative = -3.0 - PI / 2.0;
    expectedThirdDerivative = -2.0 + PI / 2.0;

    // The term cos(PI/2) introduces a small error in the function value, and in the derivative
    // values. [ cos(PI/2) is approximately 6.12323e-017 ]
    // The calculation of the function value, the first derivative and the second derivative are
    // within the numeric limits, however the error in the third derivative is larger than this
    // limit. The tolerance for the third derivative is therefore more flexible than that of the
    // other calculations.
    const double toleranceThirdDerivative = 4.0 * tolerance;

    // The exact value of the third derivative, due to the introduced error. This vaue is checked
    // separately.
    double exactThirdDerivative =
            -2.0 + PI / 2.0 - 1.0 * ( 4.0 + PI / 2.0 ) * std::cos( PI / 2.0 );

    // Initialize mathematical function.
    mission_segments::OscillatingFunctionNovak myFunctionCase2(
                boost::lambda::constant( oscillatingFunctionParameters ) );

    // Calculate function value and derivatives.
    functionValue = myFunctionCase2.evaluate( azimuthalAngle );
    firstDerivative = myFunctionCase2.computeDerivative( 1 , azimuthalAngle );
    secondDerivative = myFunctionCase2.computeDerivative( 2 , azimuthalAngle );
    thirdDerivative = myFunctionCase2.computeDerivative( 3 , azimuthalAngle );

    // Test if the computed values correspond to the expected values, within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( functionValue,
                                expectedFunctionValue,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( firstDerivative,
                                expectedFirstDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( secondDerivative,
                                expectedSecondDerivative,
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( thirdDerivative,
                                expectedThirdDerivative,
                                toleranceThirdDerivative );

    BOOST_CHECK_EQUAL( thirdDerivative, exactThirdDerivative );
}

//! Test case: request wrong derivative or integral of one of the mathematical functions.
BOOST_AUTO_TEST_CASE( wrongRequestMathematicalFunctions )
{
    // Initialize parameters of the mathematical functions.
    // Inverse polynomial parameters.
    std::pair< Eigen::Vector3d, Eigen::Vector3d > inversePolynomialParameters;
    inversePolynomialParameters.first.Zero();
    inversePolynomialParameters.second.Zero();
    double timeDepParameter = 0.0;

    // Oscillating function parameters.
    std::pair< Eigen::Vector2d, Eigen::Vector2d > OscillatingShapeParameters;
    OscillatingShapeParameters.first.Zero(  );
    OscillatingShapeParameters.second.Zero(  );

    // Initialize mathematical functions.
    mission_segments::ImprovedInversePolynomialWall myInversePolynomial(
                boost::lambda::constant( timeDepParameter ) ,
                boost::lambda::constant( inversePolynomialParameters ) );

    mission_segments::OscillatingFunctionNovak myOscillatingFunction(
                boost::lambda::constant( OscillatingShapeParameters ) );

    // Set flags.
    bool isFourthDerivativeInversePolynomial = true;
    bool isFourthDerivativeOscillatingFunction = true;
    bool isDefiniteIntegralInversePolynomial = true;
    bool isDefiniteIntegralOscillatingFunction = true;

    // Try to calculate the fourth derivative of the inverse polynomial function, which should
    // result in a runtime error.
    try
    {
        // Calculate the fourth derivative of the inverse polynomial function.
        myInversePolynomial.computeDerivative( 4 , 0.0 );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isFourthDerivativeInversePolynomial = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isFourthDerivativeInversePolynomial );

    // Try to calculate the definite integral of the inverse polynomial function, which should
    // result in a runtime error.
    try
    {
        // Calculate the definite integral of the inverse polynomial function.
        myInversePolynomial.computeDefiniteIntegral( 1 , 0.0 , 1.0 );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isDefiniteIntegralInversePolynomial = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isDefiniteIntegralInversePolynomial );

    // Try to calculate the fourth derivative of the oscillating function, which should result in a
    // runtime error.
    try
    {
        // Calculate the fourth derivative of the oscillating function.
        myOscillatingFunction.computeDerivative( 4 , 0.0 );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isFourthDerivativeOscillatingFunction = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isFourthDerivativeOscillatingFunction );

    // Try to calculate the definite integral of the oscillating function, which should result in a
    // runtime error.
    try
    {
        // Calculate the definite integral of the oscillating function.
        myOscillatingFunction.computeDefiniteIntegral( 1 , 0.0 , 1.0 );
    }

    // Catch the expected runtime error, and set the boolean flag to false.
    catch ( std::runtime_error )
    {
        isDefiniteIntegralOscillatingFunction = false;
    }

    // Check value of flag.
    BOOST_CHECK( !isDefiniteIntegralOscillatingFunction );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
