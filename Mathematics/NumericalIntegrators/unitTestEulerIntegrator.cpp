/*! \file unitTestEulerIntegrator.cpp
 *    Source file that defines the unit test for the Euler integrator included
 *    in Tudat.
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
 *    Date created      : 15 May, 2011
 *    Last modified     : 15 May, 2011
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
 *      110515    K. Kumar          File created.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include "Mathematics/NumericalIntegrators/unitTestEulerIntegrator.h"

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test implementation of Euler integrator.
bool testEulerIntegrator( )
{
    // One test.
    // Test 1: Integration of initial-value problem given on pg. 258 of
    //         (Burden and Faires, 2001).

    // Test result initialised to false.
    bool isEulerIntegratorErroneous_ = false;

    // Define tolerance for difference between computed and expected results.
    double tolerance_ = 1.0e-8;

    // Create Euler integrator.
    Euler eulerIntegrator_;

    // Set initial state in Euler integrator.
    State initialState_;
    initialState_.state.setZero( 1 );
    initialState_.state( 0 ) = 0.5;
    eulerIntegrator_.setInitialState( &initialState_ );

    // Set initial stepsize in Euler integrator.
    eulerIntegrator_.setInitialStepsize( 0.2 );

    // Set start of integration interval.
    eulerIntegrator_.setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    eulerIntegrator_.setIntegrationIntervalEnd( 2.0 );

    // Create Euler integrator test class object.
    EulerIntegratorTest eulerIntegratorTest_;

    // Set state derivative for integrator.
    eulerIntegrator_.setObjectContainingStateDerivative(
                &eulerIntegratorTest_ );

    // Execute test integration.
    eulerIntegrator_.integrate( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fabs( eulerIntegrator_.getFinalState( )->state( 0 ) - 4.8657845 )
         > tolerance_  )
    {
        isEulerIntegratorErroneous_ = true;

        std::cerr << "The computed value ( "
                  << eulerIntegrator_.getFinalState( )->state( 0 )
                  << " ) using the Euler integrator "
                  << "is not equal to the expected value: "
                  << 4.8657845 << std::endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isEulerIntegratorErroneous_;
}

//! Compute state derivative.
void EulerIntegratorTest::computeStateDerivative(
        double& time, State* pointerToState, State* pointerToStateDerivative )
{
    // Compute state derivative.
    pointerToStateDerivative->state( 0 ) = pointerToState->state( 0 )
                                           - pow( time, 2.0 ) + 1.0;
}

}

// End of file.
