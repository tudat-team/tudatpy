/*! \file unitTestEulerIntegrator.h
 *    Header file that defines the unit test for the Euler integrator included
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

#ifndef UNITTESTEULERINTEGRATOR_H
#define UNITTESTEULERINTEGRATOR_H

// Include statements.
#include <cmath>
#include <iostream>
#include "Mathematics/NumericalIntegrators/euler.h"
#include "Mathematics/NumericalIntegrators/stateDerivativeBase.h"

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test implementation of Euler integrator.
/*!
 * Tests implementation of Euler integrator.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testEulerIntegrator( );

//! Euler integrator test class.
/*!
 * Euler integrator test class.
 */
struct EulerIntegratorTest : public StateDerivativeBase
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
    void computeStateDerivative( double& time, State* pointerToState,
                                 State* pointerToStateDerivative );

protected:

private:
};

}

#endif // UNITTESTEULERINTEGRATOR_H

// End of file.
