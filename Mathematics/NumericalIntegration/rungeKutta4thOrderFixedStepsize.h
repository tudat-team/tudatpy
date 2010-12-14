/*! \file rungeKutta4thOrderFixedStepsize.h
 *    Header file that defines the 4th-order, fixed stepsize, Runge-Kutta
 *    integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegration/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dikx@student.tudelft.nl
 *
 *    Date created      : 29 September, 2010
 *    Last modified     : 29 September, 2010
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
 *      YYMMDD    author              comment
 *      100929    K. Kumar            File created.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Minor comment modifications.
 */

#ifndef RUNGEKUTTA4THORDERFIXEDSTEPSIZE_H
#define RUNGEKUTTA4THORDERFIXEDSTEPSIZE_H

// Include statements.
#include "singleStepIntegrationMethods.h"

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

    //! Perform 4th-order, fixed stepsize, Runge-Kutta integration.
    /*!
     * This function performs 4th-order, fixed stepsize, Runge-Kutta
     * integration.
     */
    void integrate( );

protected:

private:

    //! k-Coefficients
    /*!
     * k-Coefficients for 4th-order, fixed stepsize, Runge-Kutta integration
     * scheme.
     */
    VectorXd kCoefficients_[ 4 ];

    //! Modified initial state vector.
    /*!
     * Modified initial state vector, computed during a single 4th-order,
     * fixed stepsize, Runge-Kutta integration step.
     */
    VectorXd modifiedInitialStateVector_;

    //! Compute next state vector.
    /*!
     * This function computes the next state vector with respect to the state
     * vector at the start of an integration step.
     * \param stepsize Stepsize
     */
    void computeNextStateVector( double& stepsize );
};

#endif // RUNGEKUTTA4THORDERFIXEDSTEPSIZE_H

// End of file.
