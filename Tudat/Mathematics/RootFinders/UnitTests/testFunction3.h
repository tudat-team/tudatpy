/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      120318    S. Billemont      File created.
 *      120402    T. Secretin       Code-check.
 *
 *    References
 *
 *    Notes
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
