/*! \file unitTestCoefficientGenerator.h
 *    This file contains the unit test of the aerodynamic coefficient generator in Tudat.
 *
 *    Path              : /Astrodynamics/ForceModels/Aerothermodynamics/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : B.Romgens@student.tudelft.nl
 *
 *    Date created      : 4 February, 2011
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition,
 *        McGraw Hill, 2001.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd
 *        edition, AIAA Education Series, 2006
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
 *      110402    D. Dirkx          First version of file.
 *      110802    D. Dirkx          Added Apollo check.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef UNITTESTCOEFFICIENTGENERATOR_H
#define UNITTESTCOEFFICIENTGENERATOR_H

namespace unit_tests
{

//! Test coefficient Generator.
/*!
 * Tests implementation of aerodynamic coefficient generator.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testCoefficientGenerator( );

}

#endif // UNITTESTCOEFFICIENTGENERATOR_H

// End of file.
