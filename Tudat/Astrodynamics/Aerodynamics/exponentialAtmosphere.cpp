/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110224    F.M. Engelen      File created.
 *      110324    J. Melman         Added overloaded get functions.
 *      110427    F.M. Engelen      Changed input parameter to altitude, longitude and latitude.
 *      110629    F.M. Engelen      Added predefined feature.
 *      110705    F.M. Engelen      Changed to passing by reference. Changed reference values.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// The accuracy of this model could be increased by implementing different
// values for the scale height and temperature for different altitudes
// (e.g., lower, middle and upper atmosphere).
// 

#include <iostream>
#include "TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"

namespace tudat
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

        //Set specific gas constant.
        setSpecificGasConstant( physical_constants::SPECIFIC_GAS_CONSTANT_AIR );

        break;

    default:

        std::cerr << "This is not a body with a predefined exponential atmophere." << std::endl;
    }
}

} // namespace tudat
