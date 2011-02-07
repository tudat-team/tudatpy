/*! \file integratorBase.h
 *    This header file contains an abstract base class for the Integrator
 *    class included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 3
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
 *      110201    K. Kumar          First creation of code.
 *      110203    J. Melman         File checked.
 *      110207    K. Kumar          Path changed.
 */

#ifndef INTEGRATORBASE_H
#define INTEGRATORBASE_H

// Include statements.
#include "state.h"

//! An abstract base class for Integrator.
/*!
 * This is an abstract base class for the Integrator class.
 */
class IntegratorBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    IntegratorBase( ){ }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~IntegratorBase( ){ }

    //! Compute state derivative.
    /*!
     * Computes the state derivative for the numerical integrator being used.
     * \param pointerToState State given as pointer to State object.
     * \return State derivative given as pointer to State object.
     */
    virtual State* computeStateDerivative( State* pointerToState ) = 0;

protected:

private:
};

#endif // INTEGRATORBASE_H

// End of file.
