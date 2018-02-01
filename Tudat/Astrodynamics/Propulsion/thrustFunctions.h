/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_THRUSTFUNCTIONS_H
#define TUDAT_THRUSTFUNCTIONS_H

namespace tudat
{

namespace propulsion
{

//! Function to compute engine thrust from propellant mass rate and specific impulse
/*!
 * Function to compute engine thrust from propellant mass rate and specific impulse
 * \param propellantMassRate Propellant mass rate
 * \param specificImpulse Specific impulse (normalized with g0 from physical_constants namespace)
 * \return Total engine thrust
 */
double computeThrustFromSpecificImpulse(
         const double propellantMassRate, const double specificImpulse );

//! Function to compute propellant mass rate from engine thrust and specific impulse
/*!
 * Function to compute propellant mass rate from engine thrust and specific impulse
 * \param thrustMagnitude Total engine thrust
 * \param specificImpulse Specific impulse (normalized with g0 from physical_constants namespace)
 * \return Propellant mass rate
 */
double computePropellantMassRateFromSpecificImpulse(
         const double thrustMagnitude, const double specificImpulse );


} // namespace propulsion

} // namespace tudat


#endif // TUDAT_THRUSTFUNCTIONS_H
