/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      100906    S. Billemont      First creation of code.
 *      121205    K. Kumar          Migrated namespace to directory-based protocol and added
 *                                  backwards compatibility.
 *
 *    References
 *      Wolfram MathWorld, Constant:
 *          http://mathworld.wolfram.com/Constant.html (retrieved 2012/02/08).
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
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
const static double E = std::exp( 1.0 );

//! The Golden ratio \f$\approx\f$ 1.6180.
/*!
 * The golden ratio, also known as the divine proportion, golden mean, or golden section, is a
 * number often encountered when taking the ratios of distances in simple geometric figures such as
 * the pentagon, pentagram, decagon and dodecahedron.
 * \sa Wolfram MathWorld, Constant: http://mathworld.wolfram.com/GoldenRatio.html.
 */
const static double GOLDEN_RATIO = 0.5 * ( 1.0 + std::sqrt( 5.0 ) );

//! Independent root of -1, typically denoted i.
/*!
 *  Independent root of -1, typically denoted i.
 */
const static std::complex< double > COMPLEX_I = std::complex< double >( 0.0, 1.0 );

//! The constant PI \f$\approx\f$ 3.14159.
/*!
 * The constant PI, denoted \f$\pi\f$, is a real number defined as the ratio of a circle's
 * circumference C to its diameter, d = 2r.
 * \sa Wolfram MathWorld, Constant: http://mathworld.wolfram.com/Pi.html.
 */
#ifdef M_PI
const static double PI = M_PI;
#else
const static double PI = 3.141592653589793238; // 18 digits.
#endif

const static long double LONG_PI = 3.14159265358979323846264338328L;

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
ScalarType getPi( )
{
    return static_cast< ScalarType >( LONG_PI );
}


} // namespace mathematical_constants

} // namespace tudat

#endif // TUDAT_MATHEMATICAL_CONSTANTS_H
