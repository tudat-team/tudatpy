/*! \file escapePhase.h
 *    This header file contains a base class for the escape phase involved in
 *    interplanetary missions.
 *
 *    Path              : /Astrodynamics/MissionSegments/EscapeAndCapture/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 2 February, 2011
 *    Last modified     : 2 February, 2011
 *
 *    References
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
 *      110201    E. Iorfida        First creation of code.
 *
 */

#ifndef ESCAPEPHASE_H
#define ESCAPEPHASE_H

// Include statements.
#include "escapeAndCapture.h"

//! Escape phase class.
/*!
 * Escape phase class.
 */
class EscapePhase : public EscapeAndCapture
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    EscapePhase( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~EscapePhase( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param escapePhase Escape phase.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     EscapePhase& escapePhase );

protected:

private:
};

#endif // ESCAPEPHASE_H

// End of file.
