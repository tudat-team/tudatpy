/*! \file rungeKutta4thOrderFixedStepsize.h
 *    Header file that defines the 4th-order, fixed stepsize, Runge-Kutta
 *    integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 29 September, 2010
 *    Last modified     : 7 February, 2011
 *
 *    References
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
 *      100929    K. Kumar          File created.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to make use of State class.
 *      110203    J. Melman         File checked.
 *      110207    K. Kumar          Path changed; moved integrate() function
 *                                  to SingleStepIntegrationMethods.
 */

#ifndef RUNGEKUTTA4THORDERFIXEDSTEPSIZE_H
#define RUNGEKUTTA4THORDERFIXEDSTEPSIZE_H

// Include statements.
#include "singleStepIntegrationMethods.h"
#include "linearAlgebra.h"

//! 4th-order, fixed stepsize, Runge-Kutta integrator class.
/*!
 * Implementation of 4th-order, fixed stepsize, Runge-Kutta integrator.
 */
class RungeKutta4thOrderFixedStepsize : public SingleStepIntegrationMethods
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RungeKutta4thOrderFixedStepsize( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~RungeKutta4thOrderFixedStepsize( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param rungeKutta4thOrderFixedStepsize Runge-Kutta 4th-order
     *          fixed-stepsize integrator.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     RungeKutta4thOrderFixedStepsize&
                                     rungeKutta4thOrderFixedStepsize );

protected:

private:

    //! Modified current point in integration interval.
    /*!
     * Modified current point in integration interval.
     */
    double modifiedIntegrationIntervalCurrentPoint_;

    //! k-Coefficients.
    /*!
     * k-Coefficients for 4th-order, fixed stepsize, Runge-Kutta integration
     * scheme.
     */
    std::vector< VectorXd > kCoefficients_;

    //! Modified initial state.
    /*!
     * Modified initial state, computed during a single 4th-order, fixed
     * stepsize, Runge-Kutta integration step, given as a State object.
     */
    State modifiedInitialState_;

    //! Compute next state.
    /*!
     * Computes the next state with respect to the state at the start of an
     * integration step.
     * \param stepsize Stepsize.
     */
    void computeNextState_( const double& stepsize );
};

#endif // RUNGEKUTTA4THORDERFIXEDSTEPSIZE_H

// End of file.
