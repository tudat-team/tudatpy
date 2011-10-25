/*! \file orbitalElements.h
 *    This header file contains a base class for all orbital element classes in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 3
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
 *    Date created      : 20 October, 2010
 *    Last modified     : 26 October, 2010
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
 *      101020    K. Kumar          First creation of code.
 *      101022    K. Kumar          Completed constructor and set/get state functions. Added code
 *                                  comments.
 *      101026    K. Kumar          Moved set/get state functions to base class.
 */

#ifndef ORBITALELEMENTS_H
#define ORBITALELEMENTS_H

// Include statements.
#include "Astrodynamics/States/state.h"

//! Orbital elements base class.
/*!
 * Orbital elements base class.
 */
class OrbitalElements : public State
{
public:

protected:

private:
};

#endif // ORBITALELEMENTS_H

// End of file.
