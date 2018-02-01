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

#ifndef TUDAT_TEST_FUNCTION_WITH_LARGE_ROOT_DIFFERENCES_H
#define TUDAT_TEST_FUNCTION_WITH_LARGE_ROOT_DIFFERENCES_H

#include <cmath>
#include <stdexcept>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

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
