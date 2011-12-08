/*! \file unitTestRungeKuttaFehlberg78VariableStepsizeIntegrator.cpp
 *    Source file that defines the unit test for the 7(8)th-order, variable stepsize,
 *    Runge-Kutta integrator included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 22 July, 2011
 *    Last modified     : 12 September, 2011
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole,
 *          2001.
 *      Matlab RKF78 integrator, https://svn1.hosted-projects.com/cmgsoft/m_cmg
 *          /branches/nplant/ArgusTools/CILMatlab/CelestialMechMat/rkf78.m,
 *          last accessed: 26 August, 2011.
 *
 *    Notes
 *      The benchmark problem is taken from (Burden and Faires, 2001), where
 *      output data is given for a Runge-Kutta-Fehlberg 4(5)th-order
 *      integrator. Using the Matlab RKF78 integrator, the same output is
 *      generated for the given initial-value problem. The output from both is
 *      comparable.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110722    E.A.G. Heeren     File created.
 *      110803    K. Kumar          Updated unit test to use Matlab rkf78 integrator output.
 *      110811    E.A.G. Heeren     Added setTruncationErrorTolerance.
 *      110826    K. Kumar          Minor layout and comment modifications; changed
 *                                  setTruncationErrorTolerance to setRelativeErrorTolerance;
 *                                  updated benchmark data; moved includes from header file; added
 *                                  using-statement for fabs.
 *      110909    E.A.G. Heeren     Minor changes to ensure compatibility with updated
 *                                  rungeKuttaFehlberg78VariableStepsize code.
 *      110912    K. Kumar          Minor changes.
 *      111118    F.M. Engelen      Execute integration twice to check if there is no bug when
 *                                  reusing the integrator.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/NumericalIntegrators/rungeKuttaFehlberg78VariableStepsize.h"
#include "Mathematics/NumericalIntegrators/stateDerivativeBase.h"

//! Runge-Kutta-Fehlberg 7(8)th-order, variable stepsize integrator test class.
/*!
 * Runge-Kutta-Fehlberg 7(8)th-order, variable stepsize integrator test class.
 */
class RungeKuttaFehlberg78VariableStepsizeIntegratorTest : public tudat::StateDerivativeBase
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


//! Test implementation of 7(8)th-order, fixed stepsize, Runge-Kutta integrator.
int main( )
{
    // Using declarations.
    using std::fabs;
    using namespace tudat;

    // Summary of tests.
    // Test 1: Integration of initial-value problem given on pg. 278 of
    //         (Burden and Faires, 2001).

    // Test result initialised to false.
    bool isRungeKuttaFehlberg78VariableStepsizeIntegratorErroneous = false;

    // Define tolerance for difference between computed and expected results.
    double tolerance = 1.0e-14;

    // Create Runge-Kutta 7(8)th-order, fixed stepsize integrator.
    tudat::RungeKuttaFehlberg78VariableStepsize rungeKuttaFehlberg78VariableStepsizeIntegrator;

    // Set initial state in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize integrator.
    State initialState;
    initialState.state.setZero( 1 );
    initialState.state( 0 ) = 0.5;
    rungeKuttaFehlberg78VariableStepsizeIntegrator.setInitialState( &initialState );

    // Set initial stepsize in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator.setInitialStepsize( 0.01 );

    // Set truncation Error Tolerance in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator.setRelativeErrorTolerance(
                std::numeric_limits< double >::epsilon( ) );

    // Set minimum stepsize in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator.setMinimumStepsize( 1.0e-11 );

    // Set start of integration interval.
    rungeKuttaFehlberg78VariableStepsizeIntegrator.setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    rungeKuttaFehlberg78VariableStepsizeIntegrator.setIntegrationIntervalEnd( 2.0 );

    // Create Runge-Kutta-Fehlberg 7(8)th-order, variable stepsize integrator test class
    // object.
    RungeKuttaFehlberg78VariableStepsizeIntegratorTest
            rungeKuttaFehlberg78VariableStepsizeIntegratorTest;

    // Set state derivative for integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator
            .setObjectContainingStateDerivative(
                &rungeKuttaFehlberg78VariableStepsizeIntegratorTest );

    // Execute test integration.
    rungeKuttaFehlberg78VariableStepsizeIntegrator.integrate( );

    // Redo the integration.
    rungeKuttaFehlberg78VariableStepsizeIntegrator.integrate( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fabs( rungeKuttaFehlberg78VariableStepsizeIntegrator.getFinalState( )->state( 0 )
               - 5.305471950534674 ) / 5.305471950534674 > tolerance )
    {
        isRungeKuttaFehlberg78VariableStepsizeIntegratorErroneous = true;

        std::cerr << "The computed value ( "
                  << rungeKuttaFehlberg78VariableStepsizeIntegrator.getFinalState( )->state( 0 )
                  << " ) using the Runge-Kutta-Fehlberg 7(8)th-order, variable stepsize "
                  << "integrator is not equal to the expected value: "
                  <<  5.305471950534674 << std::endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isRungeKuttaFehlberg78VariableStepsizeIntegratorErroneous )
    {
        std::cerr << "testRungeKuttaFehlberg78VariableStepsizeIntegrator failed!" << std::endl;
    }

    return isRungeKuttaFehlberg78VariableStepsizeIntegratorErroneous;
}

// End of file.
