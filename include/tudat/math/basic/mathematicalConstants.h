/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wolfram MathWorld, Constant:
 *          http://mathworld.wolfram.com/Constant.html (retrieved 2012/02/08).
 */

#ifndef TUDAT_MATHEMATICAL_CONSTANTS_H
#define TUDAT_MATHEMATICAL_CONSTANTS_H

#include <complex>
#include <cmath>
#include <limits>

namespace tudat
{

namespace mathematical_constants
{

//! Constant E = exp(1) \f$\approx\f$ 2.71828.
/*!
 * The constant E is base of the natural logarithm, and is also known as Napier's constant.
 * \sa Wolfram MathWorld, Constant: http://mathworld.wolfram.com/e.html.
 */
constexpr static double E = 2.7182818284590452353602874;

//! The Golden ratio \f$\approx\f$ 1.6180.
/*!
 * The golden ratio, also known as the divine proportion, golden mean, or golden section, is a
 * number often encountered when taking the ratios of distances in simple geometric figures such as
 * the pentagon, pentagram, decagon and dodecahedron.
 * \sa Wolfram MathWorld, Constant: http://mathworld.wolfram.com/GoldenRatio.html.
 */
constexpr static double GOLDEN_RATIO = 1.6180339887498948482;

const static double EULER_MASCHERONI_CONSTANT = 0.577215664901532860606512090082;

//! Independent root of -1, typically denoted i.
/*!
 *  Independent root of -1, typically denoted i.
 */
constexpr static std::complex< double > COMPLEX_I = std::complex< double >( 0.0, 1.0 );

//! The constant PI \f$\approx\f$ 3.14159.
/*!
 * The constant PI, denoted \f$\pi\f$, is a real number defined as the ratio of a circle's
 * circumference C to its diameter, d = 2r.
 * \sa Wolfram MathWorld, Constant: http://mathworld.wolfram.com/Pi.html.
 */
#ifdef M_PI
constexpr static double PI = M_PI;
#else
constexpr static double PI = 3.14159265358979323846; // 18 digits.
#endif

constexpr static long double LONG_PI = 3.14159265358979323846264338328L;

//! Not-a-number (NaN).
/*!
 * NaN (Not a Number) is a value of the numeric data type representing an undefined or
 * unrepresentable value.
 *
 * This is a shorthand notation for std::numeric_limits<double>::signaling_NaN();
 */
#define TUDAT_NAN std::numeric_limits< double >::signaling_NaN( )

//! Function to return an integer in a floating point representation, for arbitrary
//! floating point type
/*!
 * Function to return an integer in a floating point representation, for arbitrary
 * floating point type. The function is defined as constexpr so that all operations occur
 * at compile time.
 * \param integer Integer to be represented as floating point value.
 * \return Integer in floating point representation.
 */
template< typename ScalarType  >
constexpr ScalarType getFloatingInteger( const int integer )
{
    return static_cast< ScalarType >( integer );
}

//! Function to return a rational number in a floating point representation, for arbitrary
//! floating point type
/*!
 * Function to return a rational number in a floating point representation, for arbitrary
 * floating point type. The function is defined as constexpr so that all operations occur
 * at compile time.
 * \param numerator Numerator of rational number to be represented as floating point value.
 * \param denominator Denominator of rational number to be represented as floating point value.
 * \return Integer in floating point representation.
 */
template< typename ScalarType  >
constexpr ScalarType getFloatingFraction( const int numerator, const int denominator )
{
    return static_cast< ScalarType >( numerator ) / static_cast< ScalarType >( denominator );
}

//! Function to return an the value of pi in a floating point representation, for arbitrary
//! floating point type
/*!
 * Function to return an the value of pi in a floating point representation, for arbitrary
 * floating point type. The function is defined as constexpr so that all operations occur
 * at compile time.
 * \return Pi in requested floating point representation.
 */
template< typename ScalarType  >
constexpr ScalarType getPi( )
{
    return static_cast< ScalarType >( LONG_PI );
}

} // namespace mathematical_constants

} // namespace tudat

#endif // TUDAT_MATHEMATICAL_CONSTANTS_H
