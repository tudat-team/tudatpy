/*! \file unitTestEulerIntegrator.cpp
 *    Source file that defines the unit test for the Euler integrator included
 *    in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 15 May, 2011
 *    Last modified     : 5 September, 2011
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110515    K. Kumar          File created.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Mathematics/NumericalIntegrators/euler.h"
#include "Mathematics/NumericalIntegrators/stateDerivativeBase.h"

//! Euler integrator test class.
/*!
 * Euler integrator test class.
 */
struct EulerIntegratorTest : public tudat::StateDerivativeBase
{
public:

    //! Compute state derivative.
    /*!
     * Computes state derivative. The state derivative function defined
     * corresponds to Example 1, pg. 258 in (Burden and Faires, 2001).
     * The initial-value problem is:
     * \f[
     *      y' = y - t^{ 2 } + 1
     * \f]
     * with \f$ 0 \leq t \leq 2 \f$ and \f$ y( 0 ) = 0.5 \f$.
     * \param time Time.
     * \param pointerToState Pointer to State object.
     * \param pointerToStateDerivative Computed state derivative given as a
     *          pointer to a State object.
     */
    void computeStateDerivative( double& time, tudat::State* pointerToState,
                                 tudat::State* pointerToStateDerivative )
    {
        pointerToStateDerivative->state( 0 ) = pointerToState->state( 0 )
                - std::pow( time, 2.0 ) + 1.0;
    }

protected:

private:
};

//! Test implementation of Euler integrator.
int main( )
{
    // Using declarations.
    using namespace tudat;

    // One test.
    // Test 1: Integration of initial-value problem given on pg. 258 of (Burden and Faires, 2001).

    // Test result initialised to false.
    bool isEulerIntegratorErroneous = false;

    // Define tolerance for difference between computed and expected results.
    double tolerance = 1.0e-8;

    // Create Euler integrator.
    Euler eulerIntegrator;

    // Set initial state in Euler integrator.
    State initialState_;
    initialState_.state.setZero( 1 );
    initialState_.state( 0 ) = 0.5;
    eulerIntegrator.setInitialState( &initialState_ );

    // Set initial stepsize in Euler integrator.
    eulerIntegrator.setInitialStepsize( 0.2 );

    // Set start of integration interval.
    eulerIntegrator.setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    eulerIntegrator.setIntegrationIntervalEnd( 2.0 );

    // Create Euler integrator test class object.
    EulerIntegratorTest eulerIntegratorTest;

    // Set state derivative for integrator.
    eulerIntegrator.setObjectContainingStateDerivative( &eulerIntegratorTest );

    // Execute test integration.
    eulerIntegrator.integrate( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( std::fabs( eulerIntegrator.getFinalState( )->state( 0 ) - 4.8657845 ) > tolerance )
    {
        isEulerIntegratorErroneous = true;

        std::cerr << "The computed value ( " << eulerIntegrator.getFinalState( )->state( 0 )
                  << " ) using the Euler integrator is not equal to the expected value: "
                  << 4.8657845 << std::endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isEulerIntegratorErroneous )
    {
        std::cerr << "testEulerIntegrator failed!" << std::endl;
    }

    return isEulerIntegratorErroneous;
}

// End of file.
