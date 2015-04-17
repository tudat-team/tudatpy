/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      110224    F.M. Engelen      File created.
 *      110324    J. Melman         Added overloaded get functions.
 *      110427    F.M. Engelen      Changed input parameter to altitude, longitude and latitude.
 *      110629    F.M. Engelen      Added predefined feature.
 *      110705    F.M. Engelen      Changed to passing by reference. Changed reference values.
 *
 *    References
 *
 *    Notes
 *      The accuracy of this model could be increased by implementing different values for the
 *      scale height and temperature for different altitudes (e.g., lower, middle and upper
 *      atmosphere).
 *
 */

#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"

namespace tudat
{
namespace aerodynamics
{

//! Set predefined exponential atmosphere settings.
void ExponentialAtmosphere::setPredefinedExponentialAtmosphere(
        ExponentialAtmosphere::BodiesWithPredefinedExponentialAtmospheres
        bodyWithPredefinedExponentialAtmosphere )
{
    switch( bodyWithPredefinedExponentialAtmosphere )
    {
    case earth:
        // Set local variables for Earth exponential atmosphere. Based on  lecture notes
        // Rocket Motion by Prof. Ir. B.A.C. Ambrosius, November 2009.

        // Set scale height.
        setScaleHeight( 7.200e3 );

        //Set density at zero altitude.
        setDensityAtZeroAltitude( 1.225 );

        //Set atmosphere temperature.
        setConstantTemperature( 246.0 );

        // Set specific gas constant.
        setSpecificGasConstant( physical_constants::SPECIFIC_GAS_CONSTANT_AIR );

        break;

    default:

        std::cerr << "This is not a body with a predefined exponential atmophere." << std::endl;
    }
}

} // namespace aerodynamics
} // namespace tudat
