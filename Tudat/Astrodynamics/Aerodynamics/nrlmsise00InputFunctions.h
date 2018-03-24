/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NRLMSISE00_INPUT_FUNCTIONS_H
#define TUDAT_NRLMSISE00_INPUT_FUNCTIONS_H

#include <vector>
#include <cmath>

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
 * \param altitude Altitude at which output is to be computed [m].
 * \param longitude Longitude at which output is to be computed [rad].
 * \param latitude Latitude at which output is to be computed [rad].
 * \param time Time at which output is to be computed (seconds since J2000).
 * \param solarActivityMap SolarActivityData structure
 * \param adjustSolarTime Boolean denoting whether the computed local solar time should be overidden with localSolarTime
 * input.
 * \param localSolarTime Local solar time that is used when adjustSolarTime is set to true.
 * \return NRLMSISE00Input nrlmsiseInputFunction
 */
NRLMSISE00Input nrlmsiseInputFunction( const double altitude, const double longitude,
                                       const double latitude, const double time,
                                       const tudat::input_output::solar_activity::SolarActivityDataMap& solarActivityMap,
                                       const bool adjustSolarTime = false, const double localSolarTime = 0.0 );

}  // namespace aerodynamics
}  // namespace tudat

#endif // TUDAT_NRLMSISE00_ATMOSPHERE_H_
