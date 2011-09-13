/*! \file unitTestRungeKutta4thOrderFixedStepsizeIntegrator.cpp
 *    Source file that defines the unit test for the 4th-order, fixed stepsize,
 *    Runge-Kutta integrator included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 1
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
 *    Last modified     : 17 May, 2011
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole,
 *          2001.
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
#include "Mathematics/NumericalIntegrators/unitTestRungeKutta4thOrderFixedStepsizeIntegrator.h"

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test implementation of 4th-order, fixed stepsize, Runge-Kutta integrator.
bool testRungeKutta4thOrderFixedStepsizeIntegrator( )
{
    // One test.
    // Test 1: Integration of initial-value problem given on pg. 278 of
    //         (Burden and Faires, 2001).

    // Test result initialised to false.
    bool isRungeKutta4thOrderFixedStepsizeIntegratorErroneous_ = false;

    // Define tolerance for difference between computed and expected results.
    double tolerance_ = 1.0e-8;

    // Create Runge-Kutta 4th-order, fixed stepsize integrator.
    RungeKutta4thOrderFixedStepsize rungeKutta4thOrderFixedStepsizeIntegrator_;

    // Set initial state in Runge-Kutta 4th-order, fixed stepsize integrator.
    State initialState_;
    initialState_.state.setZero( 1 );
    initialState_.state( 0 ) = 0.5;
    rungeKutta4thOrderFixedStepsizeIntegrator_
            .setInitialState( &initialState_ );

    // Set initial stepsize in Runge-Kutta 4th-order, fixed stepsize
    // integrator.
    rungeKutta4thOrderFixedStepsizeIntegrator_.setInitialStepsize( 0.2 );

    // Set start of integration interval.
    rungeKutta4thOrderFixedStepsizeIntegrator_
            .setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    rungeKutta4thOrderFixedStepsizeIntegrator_
            .setIntegrationIntervalEnd( 2.0 );

    // Create Runge-Kutta 4th-order, fixed stepsize integrator test class
    // object.
    RungeKutta4thOrderFixedStepsizeIntegratorTest
            rungeKutta4thOrderFixedStepsizeIntegratorTest_;

    // Set state derivative for integrator.
    rungeKutta4thOrderFixedStepsizeIntegrator_
            .setObjectContainingStateDerivative(
                &rungeKutta4thOrderFixedStepsizeIntegratorTest_ );

    // Execute test integration.
    rungeKutta4thOrderFixedStepsizeIntegrator_.integrate( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fabs( rungeKutta4thOrderFixedStepsizeIntegrator_
              .getFinalState( )->state( 0 ) - 5.3053630 )
         > tolerance_ )
    {
        isRungeKutta4thOrderFixedStepsizeIntegratorErroneous_ = true;

        std::cerr << "The computed value ( "
                  << rungeKutta4thOrderFixedStepsizeIntegrator_
                     .getFinalState( )->state( 0 )
                  << " ) using the Runge-Kutta 4th-order, fixed stepsize integrator "
                  << "is not equal to the expected value: "
                  << 5.3053630 << std::endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isRungeKutta4thOrderFixedStepsizeIntegratorErroneous_;
}

//! Compute state derivative.
void RungeKutta4thOrderFixedStepsizeIntegratorTest::computeStateDerivative(
        double& time, State* pointerToState, State* pointerToStateDerivative )
{
    // Compute state derivative.
    pointerToStateDerivative->state( 0 ) = pointerToState->state( 0 )
                                           - pow( time, 2.0 ) + 1.0;
}

}

// End of file.
