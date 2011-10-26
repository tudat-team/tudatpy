/*! \file state.h
 *    This header file contains a base class for all state classes in Tudat. This base class
 *    provides an interface for the rest of the Tudat library to ensure consistent implementation
 *    of state vectors, matrices etc. All state classes in Tudat are derived from this base class.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 26 October, 2010
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      The variable state is public as opposed to being protected and accessed
 *      by set and get functions to simplify the use of Eigen functions for
 *      vectors in the rest of the code, e.g. setZero(), segment().
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
 *      101026    K. Kumar          First creation of code.
 *      101110    K. Kumar          Added setState() function for matrices.
 *      101202    J. Melman         Changed path and blank line removed.
 *      110110    K. Kumar          Changed matrices to vectors.
 *      110202    K. Kumar          Removed set/get functions for state and made state a public
 *                                  variable for functionality.
 *      110204    K. Kumar          Note added about public state variable.
 *      110207    K. Kumar          Added ostream overload.
 */

#ifndef STATE_H
#define STATE_H

// Include statements.
#include <iostream>
#include "Mathematics/LinearAlgebra/linearAlgebra.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! State class.
/*!
 * Definition of State class. This is a base class for all state classes
 * defined within Tudat.
 */
class State
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~State( ) { }

    //! State.
    /*!
     * State.
     */
    VectorXd state;

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param stateObject State object.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, State& stateObject )
    { stream << "The state is set to: " << stateObject.state << std::endl; return stream; }

protected:

private:
};

}

#endif // STATE_H

// End of file.
