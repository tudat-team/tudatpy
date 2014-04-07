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

#ifndef TUDAT_TEST_FUNCTION_WITH_ZERO_ROOT_H
#define TUDAT_TEST_FUNCTION_WITH_ZERO_ROOT_H

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
 * This struct contains functions, necessary to test if the various rootfinding methods can find
 * roots that are equal to 0.
 * The test function implemented in this case is:
 * \f[
 *      f(x) = x^{2}
 * \f]
 */
struct TestFunctionWithZeroRoot : public TestFunction,
        public basic_mathematics::BasicFunction< double, double >
{
    //! Maximum order of the derivative before throwing an exception.
    unsigned int maximumDerivativeOrder;

    //! Create a function, where aMaximumDerivativeOrder is the maximum order of the derivative.
    TestFunctionWithZeroRoot( unsigned int aMaximumDerivativeOrder )
        : maximumDerivativeOrder( aMaximumDerivativeOrder )
    { }

    //! Mathematical test function.
    double evaluate( double inputValue )
    {
        // f(x) = x^2.
        return inputValue * inputValue;
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
            // Return the function value: y = inputValue^2.
            return evaluate( inputValue );
        }

        else if ( order == 1 )
        {
            // Return the first derivative function value: y = 2 * inputValue.
            return 2.0 * inputValue;
        }

        else if ( order == 2 )
        {
            // Return the second derivative function value: y = 2.
            return 2.0;
        }

        else
        {
            throw std::runtime_error(
                        "An error occured when evaluating the order of the derivative." );
        }
    }

    //! Crash on integration as root_finders should not execute these.
    double computeDefiniteIntegral( const unsigned int order, const double lowerBound,
                                    const double upperbound )
    {
        throw std::runtime_error( "The rootfinder should not evaluate integrals!" );
    }

    //! Get the expected true location of the root.
    /*!
     * Returns the expected true location of the function root, here 0.0.
     * \return True location of the root.
     */
    double getTrueRootLocation( ) { return 0.0; }
    
    //! Get the accuracy of the true location of the root.
    /*!
     * Returns the accuracy of the true location of the function root, here 1e-308.
     * \return Accuracy of the true location of the root.
     */
    double getTrueRootAccuracy( ) { return 1e-308; }

    //! Get a reasonable initial guess of the root location.
    /*!
     * Returns a reasonable initial guess for the true location of the function root, here 2.
     *
     * \return Initial guess for the true location of the function root.
     */
    double getInitialGuess( ) { return 2.0; }

    //! Get a reasonable lower boundary for the root location.
    /*!
     * Returns a reasonable lower bound for the true location of the function root, here -3.
     *
     * \return Lower bound for the true location of the function root.
     */
    double getLowerBound( ) { return -3.0; }

    //! Get a reasonable upper boundary for the root location.
    /*!
     * Returns a reasonable upper bound for the true location of the function root, here 3.
     *
     * \return Upper bound for the true location of the function root.
     */
    double getUpperBound( ) { return 3.0; }

protected:

private:
};

} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_TEST_FUNCTION_WITH_ZERO_ROOT_H
