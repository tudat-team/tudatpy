/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110627    K. Kumar          File created; ported code from previous
 *                                  predefined_gravity_field_models namespace.
 *      110701    K. Kumar          Updated references; added central gravity fields for Mercury,
 *                                  Saturn, Neptune.
 *      120327    D. Dirkx          Moved setting of predefined field to constructor.
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
 */

#include "Tudat/Astrodynamics/Gravitation/centralGravityField.h"

namespace tudat
{
namespace gravitation
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
        std::cerr << "Desired predefined central gravity field does not exist." << std::endl;
    };
}

} // namespace gravitation
} // namespace tudat
