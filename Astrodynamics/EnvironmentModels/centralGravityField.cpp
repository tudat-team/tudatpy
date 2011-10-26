/*! \file centralGravityField.cpp
 *    Source file that defines the central gravity field class in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 27 June, 2011
 *    Last modified     : 1 July, 2011
 *
 *    References
 *      Solar System Dynamics, Jet Propulsion Laboratory. Astrodynamic
 *          Constants, http://ssd.jpl.nasa.gov/?constants#ref, last accessed: 1 July, 2011.
 *      de Pater, I., Lissauer, J.J. Planetary Sciences, 2nd Edition, Cambridge
 *      University Press, Cambridge, UK, 2010.
 *      Wikipedia. Standard gravitational parameter,
 *          http://en.wikipedia.org/wiki/Standard_gravitational_parameter, last
 *          accessed: 1 July, 2011.
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
 *      110627    K. Kumar          File created; ported code from previous
 *                                  predefined_gravity_field_models namespace.
 *      110701    K. Kumar          Updated references; added central gravity fields for Mercury,
 *                                  Saturn, Neptune.
 */

// Include statements.
#include "Astrodynamics/EnvironmentModels/centralGravityField.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Tudat library namespace.
namespace tudat
{

//! Set predefined central gravity field settings.
void CentralGravityField::setPredefinedCentralGravityFieldSettings(
    BodiesWithPredefinedCentralGravityFields bodyWithPredefinedCentralGravityField )
{
    // Set degree of expansion.
    degreeOfExpansion_ = 0;

    // Set order of expansion.
    orderOfExpansion_ = 0;

    // Set reference radius to zero since for all planets, the central
    // gravity field describes a point mass representation.
    referenceRadius_ = 0.0;

    // Select body with prefined central gravity field.
    switch( bodyWithPredefinedCentralGravityField )
    {
    case sun:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 1.32712440018e20;

        break;

    case mercury:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 2.203289218e13;

        break;

    case venus:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 3.2485504415e14;

        break;

    case earth:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 3.9859383624e14;

        break;

    case moon:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 4.903686391e12;

        break;

    case mars:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.2, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 4.2828018915e13;

        break;

    case jupiter:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 1.2668579374e17;

        break;

    case saturn:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 3.793100511400001e16;

        break;

    case uranus:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 5.793943348799999e15;

        break;

    case neptune:

        // Set gravitational parameter [m^3 s^-2].
        // Reference: Mass taken from Table 1.3, pg. 6, (de Pater, 2010), value
        //            of gravitational constant taken from
        //            http://ssd.jpl.nasa.gov/?constants#ref.
        gravitationalParameter_ = 6.834733937e15;

        break;

    default:

        // Print cerr statement.
        cerr << "Desired predefined central gravity field does not exist." << endl;
    };
}

}
// End of file.
