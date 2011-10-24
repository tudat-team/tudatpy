/*! \file unitTestAerodynamicMomentAndAerodynamicForce.h
*    This header file contains the unit test for aerodynamic forces model and aerodynamic moment
*    model included in Tudat.
*
*    Path              : /Astrodynamics/ForceModels/
*    Version           : 2
*    Check status      : Checked
*
*    Checker           : F. M. Engelen
*    Affiliation       : Delft University of Technology
*    E-mail address    : F.M.Engelen@student.tudelft.nl
*
*    Checker           : D. Dirkx
*    Affiliation       : Delft University of Technology
*    E-mail address    : D.Dirkx@tudelft.nl
*
*    Date created      : 22 June, 2011
*    Last modified     : 22 August, 2011
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
*      110622    F.M. Engelen      First creation of code.
*      110822    D. Dirkx          Removed no longer necessary unit tests.
*/

#ifndef UNITTESTAERODYNAMICMOMENTANDAERODYNAMICFORCE_H
#define UNITTESTAERODYNAMICMOMENTANDAERODYNAMICFORCE_H

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test unit conversions.
/*!
 * Tests unit conversions.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testAerodynamicMomentAndAerodynamicForce( );

}

#endif // UNITTESTAERODYNAMICMOMENTANDAERODYNAMICFORCE_H

// End of file.
