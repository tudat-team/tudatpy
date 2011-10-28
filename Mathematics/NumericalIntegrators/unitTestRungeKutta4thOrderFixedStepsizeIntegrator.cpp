/*! \file unitTestRungeKutta4thOrderFixedStepsizeIntegrator.cpp
 *    Source file that defines the unit test for the 4th-order, fixed stepsize,
 *    Runge-Kutta integrator included in Tudat.
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
 *    Date created      : 17 May, 2011
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
 *      110517    K. Kumar          File created.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Mathematics/NumericalIntegrators/rungeKutta4thOrderFixedStepsize.h"
#include "Mathematics/NumericalIntegrators/stateDerivativeBase.h"

//! Runge-Kutta 4th-order, fixed stepsize integrator test class.
/*!
 *  Runge-Kutta 4th-order, fixed stepsize integrator test class.
 */
struct RungeKutta4thOrderFixedStepsizeIntegratorTest : public tudat::StateDerivativeBase
{
public:

    //! Compute state derivative.
    /*!
     * Computes state derivative. The state derivative function defined
     * corresponds to Example 3, pg. 278 in (Burden and Faires, 2001).
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

//! Test implementation of 4th-order, fixed stepsize, Runge-Kutta integrator.
int main( )
{
    // Using declarations.
    using namespace tudat;

    // One test.
    // Test 1: Integration of initial-value problem given on pg. 278 of (Burden and Faires, 2001).

    // Test result initialised to false.
    bool isRungeKutta4thOrderFixedStepsizeIntegratorErroneous = false;

    // Define tolerance for difference between computed and expected results.
    double tolerance = 1.0e-8;

    // Create Runge-Kutta 4th-order, fixed stepsize integrator.
    RungeKutta4thOrderFixedStepsize rungeKutta4thOrderFixedStepsizeIntegrator;

    // Set initial state in Runge-Kutta 4th-order, fixed stepsize integrator.
    State initialState;
    initialState.state.setZero( 1 );
    initialState.state( 0 ) = 0.5;
    rungeKutta4thOrderFixedStepsizeIntegrator.setInitialState( &initialState );

    // Set initial stepsize in Runge-Kutta 4th-order, fixed stepsize
    // integrator.
    rungeKutta4thOrderFixedStepsizeIntegrator.setInitialStepsize( 0.2 );

    // Set start of integration interval.
    rungeKutta4thOrderFixedStepsizeIntegrator.setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    rungeKutta4thOrderFixedStepsizeIntegrator.setIntegrationIntervalEnd( 2.0 );

    // Create Runge-Kutta 4th-order, fixed stepsize integrator test class
    // object.
    RungeKutta4thOrderFixedStepsizeIntegratorTest rungeKutta4thOrderFixedStepsizeIntegratorTest;

    // Set state derivative for integrator.
    rungeKutta4thOrderFixedStepsizeIntegrator.setObjectContainingStateDerivative(
                &rungeKutta4thOrderFixedStepsizeIntegratorTest );

    // Execute test integration.
    rungeKutta4thOrderFixedStepsizeIntegrator.integrate( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( std::fabs( rungeKutta4thOrderFixedStepsizeIntegrator.getFinalState( )->state( 0 )
                    - 5.3053630 ) / 5.3053630 > tolerance )
    {
        isRungeKutta4thOrderFixedStepsizeIntegratorErroneous = true;

        std::cerr << "The computed value ( "
                  << rungeKutta4thOrderFixedStepsizeIntegrator.getFinalState( )->state( 0 )
                  << " ) using the Runge-Kutta 4th-order, fixed stepsize integrator "
                  << "is not equal to the expected value: " << 5.3053630 << std::endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isRungeKutta4thOrderFixedStepsizeIntegratorErroneous )
    {
        std::cerr << "testRungeKutta4thOrderFixedStepsizeIntegrator failed!" << std::endl;
    }

    return isRungeKutta4thOrderFixedStepsizeIntegratorErroneous;
}

// End of file.
