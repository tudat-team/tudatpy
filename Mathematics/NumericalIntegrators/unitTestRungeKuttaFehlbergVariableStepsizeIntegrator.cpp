/*! \file unitTestVariableStepsizeRungeKuttaIntegrator.cpp
 *    Source file that defines the unit test for the general, variable stepsize,
 *    Runge-Kutta integrator included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Date created      : 22 July, 2011
 *    Last modified     : 10 November, 2011
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
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
 *      111021    F.M. Engelen      Added test for rkf45 and rkf56.
 *      111110    E.A.G Heeren      Minor changes.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Astrodynamics/States/state.h"
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/NumericalIntegrators/rungeKuttaFehlbergVariableStepsize.h"
#include "Mathematics/NumericalIntegrators/stateDerivativeBase.h"

//! Runge-Kutta-Fehlberg, variable stepsize integrator test class.
/*!
 * Runge-Kutta-Fehlberg, variable stepsize integrator test class.
 */
class RungeKuttaFehlbergVariableStesizeIntegratorTest : public tudat::StateDerivativeBase
{

public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RungeKuttaFehlbergVariableStesizeIntegratorTest( ) { }

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

    //! Compute the analytic answer.
    /*!
     * Analytic answer to problem, from same source as the state derivative.
     * \f[
     *      y(t) = (t+1)^{ 2 } - 0.5e^{t}
     * \f]
     * \param time Time
     * \return Analytic answer.
     */
    double computeAnalyticAnswer( double time )
    { return std::pow( ( time + 1.0 ), 2.0 ) - 0.5 * std::exp( time ); }


protected:

private:
};

//! Test implementation of the variable-stepsize, Runge-Kutta-Fehlberg integrators.
int main( )
{
    // Using declarations.
    using std::fabs;
    using tudat::mathematics::MACHINE_PRECISION_DOUBLES;
    using tudat::State;

    // Summary of tests.
    // Test 1: Integration of initial-value problem given on pg. 278 of
    //         (Burden and Faires, 2001) for the RKF78 integrator.
    // Test 2: Same as Test 1 for the RKF45 integrator.
    // Test 3: Same as Test 1 for the RKF56 integrator.

    // Test result initialized to false.
    bool isVariableStepsizeIntegratorErroneous_ = false;

    // Define tolerance for difference between computed and expected results.
    double tolerance_ = 1.0e-14;

    // Define the intergration interval end
    double intervalEnd_ = 2.0;

    // Define the initial state.
    State initialState_;
    initialState_.state.setZero( 1 );
    initialState_.state( 0 ) = 0.5;

    // Create Runge-Kutta-Fehlberg 7(8)th-order, variable stepsize integrator test class object.
    RungeKuttaFehlbergVariableStesizeIntegratorTest
            rungeKuttaFehlbergVariableStesizeIntegratorTest;

    // Test 1: RKF78 Integration of initial-value problem given on pg. 278 of
    //         (Burden and Faires, 2001).

    // Create Runge-Kutta 7(8)th-order, variable stepsize integrator.
    tudat::RungeKuttaFehlBergVariableStepsizeIntegrator
            rungeKuttaFehlberg78VariableStepsizeIntegrator_
            = tudat::RungeKuttaFehlBergVariableStepsizeIntegrator(
                    tudat::RungeKuttaFehlBergVariableStepsizeIntegrator::
                rungeKuttaFelhberg78VariableStepsize );

    // Set initial state in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.setInitialState( &initialState_ );

    // Set initial stepsize in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.setInitialStepsize( 0.01 );

    // Set truncation Error Tolerance in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.setRelativeErrorTolerance(
                MACHINE_PRECISION_DOUBLES );

    // Set minimum stepsize in Runge-Kutta Fehlberg 7(8)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.setMinimumStepsize( 1.0e-11 );

    // Set start of integration interval.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.setIntegrationIntervalEnd( intervalEnd_ );

    // Set state derivative for integrator.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.setObjectContainingStateDerivative(
                &rungeKuttaFehlbergVariableStesizeIntegratorTest );

    // Execute test integration.
    rungeKuttaFehlberg78VariableStepsizeIntegrator_.integrate( );

    double analyticResult = rungeKuttaFehlbergVariableStesizeIntegratorTest
            .computeAnalyticAnswer( intervalEnd_ );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fabs( rungeKuttaFehlberg78VariableStepsizeIntegrator_
              .getFinalState( )->state( 0 ) - analyticResult ) > tolerance_ )
    {
        isVariableStepsizeIntegratorErroneous_ = true;

        std::cerr << "The computed value ( "
                  << rungeKuttaFehlberg78VariableStepsizeIntegrator_
                     .getFinalState( )->state( 0 )
                  << " ) using the Runge-Kutta-Fehlberg 7(8)th-order, variable stepsize "
                  << "integrator is not equal to the expected value: "
                  <<  analyticResult << std::endl;
    }

    //Test 2: Test the Runge-Kutta 4(5)th-order, variable stepsize integrator.

    // Create Runge-Kutta 4(5)th-order, variable stepsize integrator.
    tudat::RungeKuttaFehlBergVariableStepsizeIntegrator
            rungeKuttaFehlberg45VariableStepsizeIntegrator_
            = tudat::RungeKuttaFehlBergVariableStepsizeIntegrator(
                    tudat::RungeKuttaFehlBergVariableStepsizeIntegrator::
                rungeKuttaFelhberg45VariableStepsize );

    // Set initial state in Runge-Kutta Fehlberg 4(5)th-order, variable stepsize integrator.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.setInitialState( &initialState_ );

    // Set initial stepsize in Runge-Kutta Fehlberg 4(5)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.setInitialStepsize( 0.01 );

    // Set truncation Error Tolerance in Runge-Kutta Fehlberg 4(5)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.setRelativeErrorTolerance(
                MACHINE_PRECISION_DOUBLES );

    // Set minimum stepsize in Runge-Kutta Fehlberg 4(5)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.setMinimumStepsize( 1.0e-11 );

    // Set start of integration interval.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.setIntegrationIntervalEnd( intervalEnd_ );

    // Set state derivative for integrator.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.setObjectContainingStateDerivative(
                &rungeKuttaFehlbergVariableStesizeIntegratorTest );

    // Execute test integration.
    rungeKuttaFehlberg45VariableStepsizeIntegrator_.integrate( );

    analyticResult = rungeKuttaFehlbergVariableStesizeIntegratorTest
            .computeAnalyticAnswer( intervalEnd_ );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fabs( rungeKuttaFehlberg45VariableStepsizeIntegrator_
              .getFinalState( )->state( 0 ) - analyticResult ) > tolerance_ )
    {
        isVariableStepsizeIntegratorErroneous_ = true;

        std::cerr << "The computed value ( "
                  << rungeKuttaFehlberg45VariableStepsizeIntegrator_
                     .getFinalState( )->state( 0 )
                  << " ) using the Runge-Kutta-Fehlberg 4(5)th-order, variable stepsize "
                  << "integrator is not equal to the expected value: "
                  <<  analyticResult << std::endl;
    }

    // Test 3: Test the Runge-Kutta 5(6)th-order, variable stepsize integrator.

    // Create Runge-Kutta 5(6)th-order, variable stepsize integrator.
    tudat::RungeKuttaFehlBergVariableStepsizeIntegrator
            rungeKuttaFehlberg56VariableStepsizeIntegrator_
            = tudat::RungeKuttaFehlBergVariableStepsizeIntegrator(
                    tudat::RungeKuttaFehlBergVariableStepsizeIntegrator::
                rungeKuttaFelhberg56VariableStepsize );

    // Set initial state in Runge-Kutta Fehlberg 5(6)th-order, variable stepsize integrator.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.setInitialState( &initialState_ );

    // Set initial stepsize in Runge-Kutta Fehlberg 5(6)th-order, variable stepsize integrator.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.setInitialStepsize( 0.01 );

    // Set truncation Error Tolerance in Runge-Kutta Fehlberg 5(6)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.setRelativeErrorTolerance(
                MACHINE_PRECISION_DOUBLES );

    // Set minimum stepsize in Runge-Kutta Fehlberg 5(6)th-order, variable stepsize
    // integrator.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.setMinimumStepsize( 1.0e-11 );

    // Set start of integration interval.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.setIntegrationIntervalStart( 0.0 );

    // Set end of integration interval.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.setIntegrationIntervalEnd( intervalEnd_ );

    // Set state derivative for integrator.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.setObjectContainingStateDerivative(
                &rungeKuttaFehlbergVariableStesizeIntegratorTest );

    // Execute test integration.
    rungeKuttaFehlberg56VariableStepsizeIntegrator_.integrate( );

    // Calculate the analytic result
    analyticResult = rungeKuttaFehlbergVariableStesizeIntegratorTest
            .computeAnalyticAnswer( intervalEnd_ );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fabs( rungeKuttaFehlberg56VariableStepsizeIntegrator_
              .getFinalState( )->state( 0 ) - analyticResult ) > tolerance_ )
    {
        isVariableStepsizeIntegratorErroneous_ = true;

        std::cerr << "The computed value ( "
                  << rungeKuttaFehlberg56VariableStepsizeIntegrator_
                     .getFinalState( )->state( 0 )
                  << " ) using the Runge-Kutta-Fehlberg 5(6)th-order, variable stepsize "
                  << "integrator is not equal to the expected value: "
                  <<  analyticResult << std::endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isVariableStepsizeIntegratorErroneous_;
}

// End of file.
