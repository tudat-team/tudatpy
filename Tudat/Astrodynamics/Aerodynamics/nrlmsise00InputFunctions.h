/*!   Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110629    F.M. Engelen      Added predefined function.
 *      110705    F.M. Engelen      Changed to passing by reference.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      151104    J. Geul           Complete rewrite
 *
 *    References
 *
 */

#ifndef TUDAT_NRLMSISE00_INPUT_FUNCTIONS_H
#define TUDAT_NRLMSISE00_INPUT_FUNCTIONS_H

#include <vector>
#include <cmath>

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h"
#include "Tudat/InputOutput/solarActivityData.h"


namespace tudat
{

namespace aerodynamics
{


//! NRLMSISE00 Input function
/*!
 * This function is used to define the input for the NRLMSISE model.
 * This function reads solar activity data and defines the input using this data.
 * \param altitude in m
 * \param longitude in RAD
 * \param latitude in RAD
 * \param time in seconds since J2000 epoch
 * \param solarActivityData SolarActivityData structure
 * \return NRLMSISE00Input nrlmsiseInputFunction
 */
NRLMSISE00Input nrlmsiseInputFunction(double altitude, double longitude,
                                              double latitude, double time,
                                              tudat::input_output::solar_activity::SolarActivityDataMap& solarActivityMap,
                                      bool adjustSolarTime = false, double localSolarTime = 0.0);

}  // namespace aerodynamics
}  // namespace tudat

#endif // TUDAT_NRLMSISE00_ATMOSPHERE_H_
