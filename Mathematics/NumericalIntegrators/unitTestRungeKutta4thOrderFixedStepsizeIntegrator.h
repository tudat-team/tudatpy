/*! \file unitTestRungeKutta4thOrderFixedStepsizeIntegrator.h
 *    Header file that defines the unit test for the 4th-order, fixed stepsize,
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
 *      110517    K. Kumar          File created.
 */

#ifndef UNITTESTRUNGEKUTTA4THORDERFIXEDSTEPSIZEINTEGRATOR_H
#define UNITTESTRUNGEKUTTA4THORDERFIXEDSTEPSIZEINTEGRATOR_H

// Include statements.
#include <cmath>
#include <iostream>
#include "rungeKutta4thOrderFixedStepsize.h"
#include "stateDerivativeBase.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test implementation of 4th-order, fixed stepsize, Runge-Kutta integrator.
/*!
 * Tests implementation of 4th-order, fixed stepsize, Runge-Kutta integrator.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testRungeKutta4thOrderFixedStepsizeIntegrator( );

//! Runge-Kutta 4th-order, fixed stepsize integrator test class.
/*!
 *  Runge-Kutta 4th-order, fixed stepsize integrator test class.
 */
class RungeKutta4thOrderFixedStepsizeIntegratorTest : public StateDerivativeBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RungeKutta4thOrderFixedStepsizeIntegratorTest( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~RungeKutta4thOrderFixedStepsizeIntegratorTest( );

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
    void computeStateDerivative( double& time, State* pointerToState,
                                 State* pointerToStateDerivative );

protected:

private:
};

}

#endif // UNITTESTRUNGEKUTTA4THORDERFIXEDSTEPSIZEINTEGRATOR_H

// End of file.
