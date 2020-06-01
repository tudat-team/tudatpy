/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_BASICELECTROMAGNETISM_H
#define TUDAT_BASICELECTROMAGNETISM_H

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{
namespace electromagnetism
{

//! Function to compute the blackbody radiation intensity at a given body temperature
/*!
 * Function to compute the blackbody radiation intensity (W/m^2) at a given body temperature.
 * \param bodyTemperature Temperature of black body emitting radiation.
 * \return Blackbodt radiation intensity at given temperature.
 */
inline double computeBlackbodyRadiationIntensity( const double bodyTemperature )
{
    return physical_constants::STEFAN_BOLTZMANN_CONSTANT *
            bodyTemperature * bodyTemperature * bodyTemperature * bodyTemperature;
}

} // namespace electromagnetism

} // namespace tudat


#endif // TUDAT_BASICELECTROMAGNETISM_H
