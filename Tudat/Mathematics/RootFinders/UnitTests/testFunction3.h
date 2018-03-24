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

#ifndef TUDAT_TEST_FUNCTION3_H
#define TUDAT_TEST_FUNCTION3_H

#include <cmath>
#include <limits>
#include <stdexcept>

#include <boost/format.hpp>             // boost::format
#include <boost/format/free_funcs.hpp>  // boost::str

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
 *      f(x) = cos(x) - x
 * \f]
 *
 * NOTE: Only derivatives up to the fifth order are defined here. You can use repetition of the
 * derivatives (see testFunction2) if you ever need higher derivatives.
 *
 */
struct TestFunction3 : public TestFunction, 
        public basic_mathematics::BasicFunction< double, double >
{
    //! Maximum order of the derivative before throwing an exception.
    unsigned int maximumDerivativeOrder;

    //! Create a function, where aMaximumDerivativeOrder is the maximum order of the derivative.
    TestFunction3( unsigned int aMaximumDerivativeOrder )
        : maximumDerivativeOrder( aMaximumDerivativeOrder )
    { }

    //! Mathematical test function.
    double evaluate( const double inputValue )
    {
        return std::cos( inputValue ) - inputValue;
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
        if ( order == 0 )
        {
            return evaluate( inputValue );
        }

        else if ( order == 1 )
        {
            return - 1.0 - std::sin( inputValue );
        }

        else if ( order == 2 )
        {
            return - std::cos( inputValue );
        }

        else if ( order == 3 )
        {
            return std::sin( inputValue );
        }

        else if ( order == 4 )
        {
            return std::cos( inputValue );
        }

        else if ( order == 5 )
        {
            return -std::sin( inputValue );
        }

        else
        {
            throw std::runtime_error( boost::str( boost::format(
                "Derivatives of order higher than 5 are not supported, requested order = %d" )
                    % order ) );
        }
        
    }

    //! Crash on integration as root_finders should not execute these.
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound,
                                    const double upperbound )
    {
        throw std::runtime_error( "The rootfinder should not evaluate integrals!" );
    }

    //! Get the expected true location of the root.
    double getTrueRootLocation( )
    {
        // Result computed using arbitrary precision in Mathematica, way more precise than IEEE
        // double: FindRoot[Cos[x]-x, {x, 3}, WorkingPrecision -> 48];
        return 0.739085133215160641655312087673873404013411758901;
    }
    
    //! Get the accuracy of the true location of the root.
    double getTrueRootAccuracy( )
    {
        return getTrueRootLocation( ) * std::numeric_limits< double >::epsilon( );
    }

    //! Get a reasonable initial guess of the root location = -2.
    double getInitialGuess( ) { return -2.0; }

    //! Get a reasonable lower boundary for the root location = -1.
    double getLowerBound( ) { return -1.0; }

    //! Get a reasonable upper boundary for the root location = 2.
    double getUpperBound( ) { return 2.0; }

protected:

private:
};

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_FUNCTION3_H
