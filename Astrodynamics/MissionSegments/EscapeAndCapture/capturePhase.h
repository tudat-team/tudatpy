/*! \file capturePhase.h
 *    This header file contains a base class for the capture phase involved in
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

#ifndef CAPTUREPHASE_H
#define CAPTUREPHASE_H

// Include statements.
#include "escapeAndCapture.h"

//! Capture phase class.
/*!
 * Capture phase class.
 */
class CapturePhase : public EscapeAndCapture
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CapturePhase( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CapturePhase( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param capturePhase Capture phase.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     CapturePhase& capturePhase );

protected:

private:
};

#endif // CAPTUREPHASE_H

// End of file.
