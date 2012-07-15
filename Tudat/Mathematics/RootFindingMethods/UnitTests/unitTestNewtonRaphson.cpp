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
 *      110111    E. Iorfida        First creation of the code.
 *                                  The code is tested with the function: f(x)= x^2 - 3.
 *      110111    K. Kumar          Updated to use address of global  functions instead of
 *                                  pointers for set functions; aligned code as required for
 *                                  namespaces; minor comment changes.
 *      110119    K. Kumar          Updated code to work with adaptor and abstract base
 *                                  implementation so that pointer-to-member functions are not
 *                                  required; filename changed; added cerr statements.
 *      110120    E. Iorfida        Added necessary class that contains functions, and related
 *                                  code, to allow a directly test with adaptor.
 *      110120    K. Kumar          Added global functions test; updated comments; modified layout.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120712    P. Musegaas       Changed absolute tolerance into a safe variant of relative
 *                                  tolerance. Added new unit tests for this.
 *
 *    References
 *
 *
 *    Notes
 *      The current implementation of the tolerance (a relative tolerance) may not be ideal for all
 *      applications (i.e. not good for values close to 0.0). It was selected because the old
 *      implementation was not suited for functions whose root might differ orders of magnitude.
 *      This version is safe for all applications though. Many iterations may be required if one
 *      searches for roots close to zero but not typically equal to zero.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
namespace unit_tests
{

//! Struct for NewtonRaphson unit test code.
/*!
 * This struct contains functions, necessary to test NewtonRaphson method.
 */
struct NewtonRaphsonTest
{
public:

    //! Mathematical test function.
    /*!
     * Mathematical test function used by the Newton-Raphson algorithm.
     * \param inputValue Input value.
     */
    double computeTestFunction( double& inputValue ) { return inputValue * inputValue - 3.0; }

    //! First-derivative of mathematical test function.
    /*!
     * First-derivative of mathematical test function used by the
     * Newton-Raphson algorithm.
     * \param inputValue Input value.
     */
    double computeFirstDerivativeTestFunction( double& inputValue ) { return 2.0 * inputValue; }

protected:

private:
};

//! Complute global mathematical test function.
/*!
 * Computes global test function:
 * \f[
 *      y = x^{2} - 3
 * \f]
 * This function is used to illustrate the use of global free functions with the Newton-Raphon
 * root-finder.
 * \param inputValue Input value (\f$x\f$).
 * \return Computed value (\f$y\f$).
 * \sa computeGlobalFirstDerivativeTestFunction().
 */
double computeGlobalTestFunction( double& inputValue )
{
    return inputValue * inputValue - 3.0;
}

//! Compute global first-derivative mathematical test function.
/*!
 * Computes global first-derivative of test function:
 * \f[
 *      \frac{dy}{dx} = y' = 2 * x
 * \f]
 * This function is used to illustrate the use of global free functions with the Newton-Raphon
 * root-finder.
 * \param inputValue Input value (\f$x\f$).
 * \return Computed value (\f$y'\f$).
 * \sa computeGlobalTestFunction().
 */
double computeGlobalFirstDerivativeTestFunction( double& inputValue ) { return 2.0 * inputValue; }

//! Compute zero-root function.
/*!
 * A simple function whose root is zero.
 * \param inputValue Input value.
 * \return Computed value.
 */
double computeZeroRootFunction( double& inputValue ) { return inputValue * inputValue; }

//! Compute first-derivative of zero-root function.
/*!
 * Computes first-derivative of the simple function with zero root.
 */
double computeFirstDerivativeZeroRootFunction( double& inputValue ) { return 2.0 * inputValue; }

//! Struct for testing a function with large differences in roots.
/*!
 * This struct contains functions to test if the root finder converges correctly to a function
 * whose roots vary a lot. Similar to the eccentricity finding functions in a gravity assist.
 */
struct NewtonRaphsonLargeRootDifferencesTest
{
public:

    //! Constructor with immediate definition of parameters.
    /*!
     * Constructor that sets all the parameters in the eccentricity finding functions for use in the
     * Newton-Raphson root-finder.
     */
    NewtonRaphsonLargeRootDifferencesTest ( const double incomingSemiMajorAxis,
                                            const double outgoingSemiMajorAxis,
                                            const double bendingAngle )
        : incomingSemiMajorAxis_( incomingSemiMajorAxis),
          outgoingSemiMajorAxis_( outgoingSemiMajorAxis ),
          bendingAngle_ ( bendingAngle )
    { }

    //! Compute incoming eccentricity function.
    /*!
     * Computes incoming eccentricity function. This function is used by the Newton-Raphson root-
     * finder to find the incoming eccentricity that matches the bending angle required in the
     * gravity assist.
     * \param incomingEccentricity Incoming eccentricity.
     * \return Incoming eccentricity root finding function value.
     * \sa NewtonRaphson().
     */
    double computeIncomingEccentricityFunction( double& incomingEccentricity )
    {
        return std::asin( 1.0 / incomingEccentricity )
                + std::asin( 1.0 / ( 1.0 - incomingSemiMajorAxis_ / outgoingSemiMajorAxis_ *
                                     ( 1.0 - incomingEccentricity ) ) ) - bendingAngle_;
    }

    //! Compute first-derivative of the incoming eccentricity function.
    /*!
     * Computes the first-derivative of the incoming eccentricity function. This function is used
     * by the Newton-Raphson root-finder to find the incoming eccentricity that matches the bending
     * angle required in the gravity assist.
     * \param incomingEccentricity Incoming eccentricity.
     * \return Incoming eccentricity root finding function first-derivative value.
     * \sa NewtonRapshon().
     */
    double computeFirstDerivativeIncomingEccentricityFunction( double& incomingEccentricity )
    {
        const double eccentricitySquareMinusOne_ =
                incomingEccentricity * incomingEccentricity - 1.0;
        const double semiMajorAxisRatio_ = incomingSemiMajorAxis_ / outgoingSemiMajorAxis_ ;
        const double bParameter_ = 1.0 - semiMajorAxisRatio_ * ( 1.0 - incomingEccentricity );

        return -1.0 / ( incomingEccentricity * std::sqrt( eccentricitySquareMinusOne_ ) ) -
               semiMajorAxisRatio_ / ( bParameter_ * std::sqrt( bParameter_ * bParameter_ - 1.0 ) );
    }

protected:

private:

    //! Semi-major axis of the incoming hyperbolic leg
    const double incomingSemiMajorAxis_;

    //! Semi-major axis of the outgoing hyperbolic leg.
    const double outgoingSemiMajorAxis_;

    //! Bending angle between the excess velocities.
    const double bendingAngle_;
};

BOOST_AUTO_TEST_SUITE( test_newton_raphson )

//! Test if Newton-Raphson root-finder works correctly using global functions.
BOOST_AUTO_TEST_CASE( testNewtonRaphsonWithGlobalFunctions )
{
    // Set expected root.
    const double expectedRoot = std::sqrt( 3.0 );

    // Declare new Newton-Raphson object.
    boost::shared_ptr< tudat::NewtonRaphson > newtonRaphson
            = boost::make_shared< tudat::NewtonRaphson >( );

    // Set values for the implementation of the code.
    newtonRaphson->setRelativeTolerance( 1.0e-15 );
    newtonRaphson->setInitialGuessOfRoot( 5.0 );

    // Set mathematical functions.
    newtonRaphson->setMathematicalFunction( &computeGlobalTestFunction );
    newtonRaphson->setFirstDerivativeMathematicalFunction(
                &computeGlobalFirstDerivativeTestFunction );

    // Compute root.
    newtonRaphson->execute( );

    // Check if computed root matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedRoot, newtonRaphson->getComputedRootOfFunction( ),
                                newtonRaphson->getRelativeTolerance( ) );

}

//! Test if Newton-Raphson root-finder works correctly using member functions.
BOOST_AUTO_TEST_CASE( testNewtonRaphsonWithMemberFunctions )
{
    // Set expected root.
    const double expectedRoot = std::sqrt( 3.0 );

    // Declare new Newton-Raphson object.
    boost::shared_ptr< tudat::NewtonRaphson > newtonRaphson
            = boost::make_shared< tudat::NewtonRaphson >( );

    // Set values for the implementation of the code.
    newtonRaphson->setRelativeTolerance( 1.0e-15 );
    newtonRaphson->setInitialGuessOfRoot( 5.0 );

    // Declare NewtonRaphsonAdaptor object.
    tudat::NewtonRaphsonAdaptor< NewtonRaphsonTest > newtonRaphsonAdaptor_;

    // Set adaptor class object and member functions.
    newtonRaphson->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptor_ );
    newtonRaphsonAdaptor_.setPointerToFunction(
                &NewtonRaphsonTest::computeTestFunction );
    newtonRaphsonAdaptor_.setPointerToFirstDerivativeFunction(
                &NewtonRaphsonTest::computeFirstDerivativeTestFunction );

    // Compute root.
    newtonRaphson->execute( );

    // Check if computed root matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedRoot, newtonRaphson->getComputedRootOfFunction( ),
                                newtonRaphson->getRelativeTolerance( ) );

}

//! Test if Newton-Raphson root-finder works correctly with functions whose root is zero.
BOOST_AUTO_TEST_CASE( testNewtonRaphsonWithZeroRoot )
{
    // Declare tolerance.
    const double tolerance = std::numeric_limits< double >::min( );

    // Declare new Newton-Raphson object.
    boost::shared_ptr< tudat::NewtonRaphson > newtonRaphson
            = boost::make_shared< tudat::NewtonRaphson >( );

    // Set values for the implementation of the code.
    newtonRaphson->setRelativeTolerance( 1.0e-15 );
    newtonRaphson->setZeroRepresentation( 1.0e-20 );
    newtonRaphson->setInitialGuessOfRoot( 5.0 );

    // Set adaptor class object and member functions.
    newtonRaphson->setMathematicalFunction( &computeZeroRootFunction );
    newtonRaphson->setFirstDerivativeMathematicalFunction(
                &computeGlobalFirstDerivativeTestFunction );

    // Compute root.
    newtonRaphson->execute( );

    // Check if computed root matches expected value.
    BOOST_CHECK_SMALL( newtonRaphson->getComputedRootOfFunction( ), tolerance );
}

//! Test if Newton-Raphson root-finder works correctly with functions having very different roots.
BOOST_AUTO_TEST_CASE( testNewtonRaphsonWithLargeDifferencesInRoots )
{
    // Declare tolerance.
    const double tolerance = 1.0e-10;

    // Declare expected roots.
    double expectedRootLow = 1.00000000793634;
    double expectedRootHigh = 7937.3386333591;

    // Very similar semi major axes and bending angles have arisen in trajectory optimization of
    // Cassini, hence although they are not realistic (especially within the patched conics
    // framework), they need to be calculated correctly. In general much more constraining
    // situations can be thought of. This unit test can be upgraded.
    double incomingSemiMajorAxisLow = -3.24859999867635e18;
    double outgoingSemiMajorAxisLow = -3248600.0;
    double bendingAngleLow = 1.5707963267949;

    double incomingSemiMajorAxisHigh = -3248600.0;
    double outgoingSemiMajorAxisHigh = -3.24859999867635e18;
    double bendingAngleHigh = 1.5707963267949;

    // Instantiate the low- and high-case classes.
    NewtonRaphsonLargeRootDifferencesTest functionsLow( incomingSemiMajorAxisLow,
                                                        outgoingSemiMajorAxisLow,
                                                        bendingAngleLow );
    NewtonRaphsonLargeRootDifferencesTest functionsHigh( incomingSemiMajorAxisHigh,
                                                         outgoingSemiMajorAxisHigh,
                                                         bendingAngleHigh );

    // Declare new Newton-Raphson objects.
    boost::shared_ptr< tudat::NewtonRaphson > newtonRaphsonLow
            = boost::make_shared< tudat::NewtonRaphson >( );

    boost::shared_ptr< tudat::NewtonRaphson > newtonRaphsonHigh
            = boost::make_shared< tudat::NewtonRaphson >( );

    // Set values for the implementation of the code.
    newtonRaphsonLow->setRelativeTolerance( tolerance );
    newtonRaphsonLow->setInitialGuessOfRoot( 1.0 + 1.0e-10 );
    newtonRaphsonHigh->setRelativeTolerance( tolerance );
    newtonRaphsonHigh->setInitialGuessOfRoot( 1.0 + 1.0e-2 );

    // Declare NewtonRaphsonAdaptor objects.
    tudat::NewtonRaphsonAdaptor< NewtonRaphsonLargeRootDifferencesTest > newtonRaphsonAdaptorLow;
    newtonRaphsonAdaptorLow.setClass( &functionsLow );
    tudat::NewtonRaphsonAdaptor< NewtonRaphsonLargeRootDifferencesTest > newtonRaphsonAdaptorHigh;
    newtonRaphsonAdaptorHigh.setClass( &functionsHigh );

    // Set adaptor class object and member functions.
    newtonRaphsonLow->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptorLow );
    newtonRaphsonAdaptorLow.setPointerToFunction(
                &NewtonRaphsonLargeRootDifferencesTest::computeIncomingEccentricityFunction );
    newtonRaphsonAdaptorLow.setPointerToFirstDerivativeFunction(
                &NewtonRaphsonLargeRootDifferencesTest::
                        computeFirstDerivativeIncomingEccentricityFunction );
    newtonRaphsonHigh->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptorHigh );
    newtonRaphsonAdaptorHigh.setPointerToFunction(
                &NewtonRaphsonLargeRootDifferencesTest::computeIncomingEccentricityFunction );
    newtonRaphsonAdaptorHigh.setPointerToFirstDerivativeFunction(
                &NewtonRaphsonLargeRootDifferencesTest::
                        computeFirstDerivativeIncomingEccentricityFunction );

    // Compute root.
    newtonRaphsonLow->execute( );
    newtonRaphsonHigh->execute( );

    // Check if computed root matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedRootLow, newtonRaphsonLow->getComputedRootOfFunction( ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedRootHigh, newtonRaphsonHigh->getComputedRootOfFunction( ),
                                tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
