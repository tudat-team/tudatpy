/*! \file state.h
 *    This header file contains a base class for all state classes in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
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
 *    Date created      : 26 Checked, 2010
 *    Last modified     : 02 December, 2010
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
 *      YYMMDD    author        comment
 *      101026    K. Kumar      First creation of code.
 *      101202    J. Melman     Changed path and blank line removed.
 */

#ifndef STATE_H
#define STATE_H

// Include statements.
#include "linearAlgebra.h"

//! State class.
/*!
 * Definition of State class.
 */
class State
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    State( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~State( );

    //! Set state.
    /*!
     * This function sets the state.
     * \param state State given as a vector.
     */
    void setState( VectorXd& state );

    //! Get state.
    /*!
     * This function returns the state.
     * \return State.
     */
    VectorXd& getState( );

protected:

    //! State.
    /*!
     * State.
     */
    VectorXd state_;

private:
};

#endif // STATE_H

// End of file.
