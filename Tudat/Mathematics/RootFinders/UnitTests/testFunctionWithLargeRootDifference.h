/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      120719    P. Musegaas       File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_TEST_FUNCTION_WITH_LARGE_ROOT_DIFFERENCES_H
#define TUDAT_TEST_FUNCTION_WITH_LARGE_ROOT_DIFFERENCES_H

#include <cmath>
#include <stdexcept>

#include "Tudat/Mathematics/RootFinders/UnitTests/testFunction.h"
#include "Tudat/Mathematics/BasicMathematics/basicFunction.h"

namespace tudat
{
namespace unit_tests
{

//! Test function for the root-finders.
/*!
 * This struct contains a function which potentially has a big difference in the root. Based on the
 * eccentricity finding function in the gravity assist method.
 */
struct TestFunctionWithLargeRootDifference : public TestFunction,
        public basic_mathematics::BasicFunction< double, double >
{
    //! Maximum order of the derivative before throwing an exception.
    unsigned int maximumDerivativeOrder;

    //! Create a function, where aMaximumDerivativeOrder is the maximum order of the derivative.
    TestFunctionWithLargeRootDifference( const unsigned int aMaximumDerivativeOrder,
                                         const double anIncomingSemiMajorAxis,
                                         const double anOutgoingSemiMajorAxis,
                                         const double aBendingAngle )
        : maximumDerivativeOrder( aMaximumDerivativeOrder ),
          incomingSemiMajorAxis( anIncomingSemiMajorAxis ),
          outgoingSemiMajorAxis( anOutgoingSemiMajorAxis ),
          bendingAngle( aBendingAngle )
    { }

    //! Mathematical test function.
    double evaluate( const double inputValue )
    {
        return std::asin( 1.0 / inputValue )
                + std::asin( 1.0 / ( 1.0 - incomingSemiMajorAxis / outgoingSemiMajorAxis *
                                     ( 1.0 - inputValue ) ) ) - bendingAngle;
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
            // Return the function value.
            return evaluate( inputValue );
        }

        else if ( order == 1 )
        {
            // Return the first derivative function value.
            const double eccentricitySquareMinusOne = inputValue * inputValue - 1.0;
            const double semiMajorAxisRatio = incomingSemiMajorAxis / outgoingSemiMajorAxis;
            const double bParameter = 1.0 - semiMajorAxisRatio * ( 1.0 - inputValue );

            return -1.0 / ( inputValue * std::sqrt( eccentricitySquareMinusOne ) ) -
                semiMajorAxisRatio / ( bParameter * std::sqrt( bParameter * bParameter - 1.0 ) );
        }

        else
        {
            throw std::runtime_error(
                        "An error occured when evaluating the order of the derivative." );
        }
    }

    //! Crash on integration as root_finders should not execute these.
    double computeDefiniteIntegral( unsigned int order, double lowerBound, double upperbound )
    {
        throw std::runtime_error( "The rootfinder should not evaluate integrals!" );
    }

    //! Get the expected true location of the root.
    /*!
     * Not implemented.
     */
    double getTrueRootLocation( ) { return TUDAT_NAN; }
    
    //! Get the accuracy of the true location of the root.
    /*!
     * Not implemented.
     */
    double getTrueRootAccuracy( ) { return TUDAT_NAN; }

    //! Get a reasonable initial guess of the root location.
    /*!
     * Not implemented.
     */
    double getInitialGuess( ) { return TUDAT_NAN; }

    //! Get a reasonable lower boundary for the root location.
    /*!
     * Not implemented.
     */
    double getLowerBound( ) { return TUDAT_NAN; }

    //! Get a reasonable upper boundary for the root location.
    /*!
     * Not implemented.
     */
    double getUpperBound( ) { return TUDAT_NAN; }

    //! Semi-major axis of the incoming hyperbolic leg
    const double incomingSemiMajorAxis;

    //! Semi-major axis of the outgoing hyperbolic leg.
    const double outgoingSemiMajorAxis;

    //! Bending angle between the excess velocities.
    const double bendingAngle;

protected:

private:
};

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_FUNCTION_WITH_LARGE_ROOT_DIFFERENCES_H
