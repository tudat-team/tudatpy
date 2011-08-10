/*! \file stateDerivativeBase.h
 *    This header file contains an abstract base class to compute state
 *    derivatives, for use with the numerical integrators included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 1 February, 2011
 *    Last modified     : 10 August, 2011
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
 *      110201    K. Kumar          First creation of code.
 *      110203    J. Melman         File checked.
 *      110516    K. Kumar          Renamed file and class.
 *      110810    J. Leloux         Corrected doxygen documentation, deleted
 *                                  unused variable.
 */

#ifndef STATEDERIVATIVEBASE_H
#define STATEDERIVATIVEBASE_H

// Include statements.
#include "state.h"

//! An abstract base class for state derivative functions.
/*!
 * An abstract base class for state derivative functions.
 */
class StateDerivativeBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    StateDerivativeBase( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~StateDerivativeBase( );

    //! Compute state derivative.
    /*!
     * Computes the state derivative for the numerical integrator being used.
     * \param independentVariable Value of independent variable.
     * \param pointerToState State given as pointer to State object.
     * \param pointerToStateDerivative State derivative given as pointer to
     *          State object. The computed  state derivative is stored in the
     *          object that this pointer points to.
     */
    virtual void computeStateDerivative( double& independentVariable,
                                         State* pointerToState,
                                         State* pointerToStateDerivative ) = 0;

protected:

private:
};

#endif // STATEDERIVATIVEBASE_H

// End of file.
