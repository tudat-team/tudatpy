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

#ifndef TUDAT_TEST_FUNCTION2_H
#define TUDAT_TEST_FUNCTION2_H

#include <cmath>
#include <limits>
#include <stdexcept>

#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction.h"
#include "Tudat/Mathematics/BasicMathematics/basicFunction.h"

namespace tudat
{
namespace unit_tests
{

//! Test function for the root_finders.
/*!
 * This struct contains functions, necessary to test the various rootfinding methods.
 *
 * The test function implemented in this case is:
 * \f[
 *      f(x) = sin(x)
 * \f]
 */
struct TestFunction2 : public TestFunction,
        public basic_mathematics::BasicFunction< double, double >
{
    //! Maximum order of the derivative before throwing an exception.
    unsigned int maximumDerivativeOrder;

    //! Create a function, where aMaximumDerivativeOrder is the maximum order of the derivative.
    TestFunction2( const unsigned int aMaximumDerivativeOrder )
        : maximumDerivativeOrder( aMaximumDerivativeOrder )
    { }

    //! Mathematical test function.
    double evaluate( const double inputValue )
    {
        return std::sin( inputValue );
    }

    //! Derivatives of mathematical test function.
    double computeDerivative( const unsigned int order, const double inputValue )
    { 
        // Sanity check.
        if ( order > maximumDerivativeOrder )
        {
            throw std::runtime_error( "The rootfinder should not evaluate higher derivatives!" );
        }

        // Return the analytical expression for the derivatives.
        if ( order % 4 == 0 )
        {
            return evaluate( inputValue );
        }

        else if ( order % 4 == 1 )
        {
            return std::cos( inputValue );
        }

        else if ( order % 4 == 2 )
        {
            return - std::sin( inputValue );
        }

        else // if ( order % 4 == 3 ) The only other option.
        {
            return - std::cos( inputValue );
        }
    }

    //! Crash on integration as root_finders should not execute these.
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound,
                                    const double upperbound )
    {
        throw std::runtime_error( "The rootfinder should not evaluate integrals!" );
    }

    //! Get the expected true location of the root = sqrt(3).
    double getTrueRootLocation( )
    {
        // Result computed using arbitrary precision in Mathematica, much more precise than IEEE
        // double: FindRoot[Sin[x], {x, 3}, WorkingPrecision -> 48];
        return 3.14159265358979323846264338327950288419716939938;
    }
    
    //! Get the accuracy of the true location of the root.
    double getTrueRootAccuracy( )
    {
        return getTrueRootLocation( ) * std::numeric_limits< double >::epsilon( );
    }

    //! Get a reasonable initial guess of the root location = 3.
    double getInitialGuess( ) { return 3.0;  }

    //! Get a reasonable lower boundary for the root location = 2.
    double getLowerBound( ) { return 2.0; }

    //! Get a reasonable upper boundary for the root location = 4.
    double getUpperBound( ) { return 4.0; }

protected:

private:
};

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_FUNCTION2_H
