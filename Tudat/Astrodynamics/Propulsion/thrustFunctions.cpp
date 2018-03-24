/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "Tudat/Astrodynamics/Propulsion/thrustFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

namespace propulsion
{

//! Function to compute engine thrust from propellant mass rate and specific impulse
double computeThrustFromSpecificImpulse(
         const double propellantMassRate, const double specificImpulse )
{
    return  physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * specificImpulse * propellantMassRate;
}

//! Function to compute propellant mass rate from engine thrust and specific impulse
double computePropellantMassRateFromSpecificImpulse(
         const double thrustMagnitude, const double specificImpulse )
{
    return  thrustMagnitude / ( physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION * specificImpulse );
}


} // namespace propulsion

} // namespace tudat

