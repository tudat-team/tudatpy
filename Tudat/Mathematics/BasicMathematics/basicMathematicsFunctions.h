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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Spiegel, M.R., Stephens, L.J. Statistics, Fourth Edition, Schaum's Outlines, McGraw-Hill,
 *          2008.
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
 */

#ifndef TUDAT_BASIC_MATHEMATICS_FUNCTIONS_H
#define TUDAT_BASIC_MATHEMATICS_FUNCTIONS_H

#include <boost/random/mersenne_twister.hpp>

namespace tudat
{
namespace basic_mathematics
{

//! Random number generator typedef.
/*!
 * Random number generator typedef. This can be modified to any other Boost random number
 * generator type (http://www.boost.org/doc/libs/1_47_0/doc/html/boost_random/reference.html).
 */
typedef boost::mt19937 GlobalRandomNumberGeneratorType;

//! Get global random number generator.
/*!
 * Returns global random number generator. The default seed is set to the current time.
 * \return Global random number generator.
 */
GlobalRandomNumberGeneratorType& getGlobalRandomNumberGenerator( );

//! Compute modulo of floating-point number (default double).
/*!
 * Computes the remainder of division of one floating-point number by another. The modulo
 * computation is based on the mathematical definition of congruence, which is different from the
 * implementation of std::fmod() in the cmath standard library. For a description of congruence
 * see: http://mathworld.wolfram.com/Congruence.html.
 * The remainder is in the range [ 0, divisor ).
 * \param dividend Number to be divided.
 * \param divisor Number that is divided by.
 * \return Remainder of division of dividend by divisor.
 */
template< typename ScalarType = double >
ScalarType computeModulo( const ScalarType dividend, const ScalarType divisor )
{
    return dividend - divisor * std::floor( dividend / divisor );
}

//! Compute modulo of floating-point number (default double).
/*!
 * Computes the remainder of division of one floating-point number by another. The modulo
 * computation is based on the mathematical definition of congruence, which is different from the
 * implementation of std::fmod() in the cmath standard library. For a description of congruence
 * see: http://mathworld.wolfram.com/Congruence.html. The remainder is in the range [ 0, divisor ).
 * This function also returns (by reference) the number of times divisor goes into dividend, i.e. the division from which
 * the moduloValue is the remainder.
 * \param dividend Number to be divided.
 * \param divisor Number that is divided by.
 * \param moduloValue Remainder of division of dividend by divisor (returned by reference).
 * \param numberOfDivisors Number of times divisor goes into dividend (returned by reference).
 */
template< typename ScalarType >
inline void computeModuloAndRemainder( const ScalarType dividend, const ScalarType divisor,
                                       ScalarType& moduloValue, int& numberOfDivisors )
{
    numberOfDivisors = std::floor( dividend / divisor );
    moduloValue = dividend - divisor * static_cast< ScalarType >( numberOfDivisors );
}


//! Raise floating point variable to integer power.
template< typename ScalarType >
ScalarType raiseToIntegerPower( const ScalarType baseValue,
                            const int integerPower )
{
    // Declare local variable.
    // Declare result of raising base to integer power.
    // Initialise with value.
    ScalarType resultOfRaisingBaseToIntegerPower = 1;
    // Declare absolute value of integerPower.
    int absoluteValueOfIntegerPower
            = std::abs( integerPower );
    // Declare copy of base value.
    ScalarType copyOfBaseValue = baseValue;

    // Compute the result here using exponentiation by squares.
    // Stop loop when absolute value of integer power is equal to zero.
    while ( absoluteValueOfIntegerPower )
    {
        // Check that absolute value of integer power.
        if ( absoluteValueOfIntegerPower & 1 )
        {
            // Compute intermediate result.
            resultOfRaisingBaseToIntegerPower *= copyOfBaseValue;
        }

        // Divide integer power by two.
        absoluteValueOfIntegerPower >>= 1;

        // Square base value.
        copyOfBaseValue *= copyOfBaseValue;
    }

    // Check if sign of integerPower is negative.
    if ( integerPower < 0 )
    {
        // Switch sign of result.
        resultOfRaisingBaseToIntegerPower
                = 1.0 / resultOfRaisingBaseToIntegerPower;
    }

    // Return result of raising base to integer power.
    return resultOfRaisingBaseToIntegerPower;
}

} // namespace basic_mathematics
} // namespace tudat


#endif // TUDAT_BASIC_MATHEMATICS_FUNCTIONS_H
