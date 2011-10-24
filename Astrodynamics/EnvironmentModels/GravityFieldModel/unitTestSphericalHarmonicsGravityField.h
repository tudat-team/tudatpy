/*! \file unitTestSphericalHarmonicsGravityField.h
 *   Header file for a unit test that tests the implementation of the spherical
 *   harmonics gravity field class in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/GravityFieldModel/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 7 January, 2011
 *    Last modified     : 13 January, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
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
 *      110107    K. Kumar          First creation of code.
 *      110113    K. Kumar          Minor updates to match protocol.
 */

#ifndef UNITTESTSPHERICALHARMONICSGRAVITYFIELD_H
#define UNITTESTSPHERICALHARMONICSGRAVITYFIELD_H

//! Namespace for all unit tests.
/*!
 * Namespace containing all unit tests.
 */
namespace unit_tests
{

//! Test of implementation of spherical harmonics gravity field class.
/*!
 * Test of implementation of spherical harmonics gravity field class.
 * \return Boolean indicating success of test
 * ( false = successful; true = failed ).
 */
bool testSphericalHarmonicsGravityField( );

}

#endif // UNITTESTSPHERICALHARMONICSGRAVITYFIELD_H

// End of file.
