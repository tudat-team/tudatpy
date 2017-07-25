/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_PROPAGATOR_H
#define TUDAT_JSONINTERFACE_PROPAGATOR_H

#include <Tudat/SimulationSetup/PropagationSetup/propagationSettings.h>

#include "Tudat/External/JsonInterface/jsonInterface.h"

namespace tudat
{

namespace propagators
{

//! Map of `IntegratedStateType`s supported by `json_interface`.
static std::map< std::string, IntegratedStateType > integratedStateTypes =
{
    { "hybrid",               hybrid },
    { "translational",	      translational_state },
    { "rotational",           rotational_state },
    { "bodyMass",             body_mass_state },
    { "relativisticTimeRate", relativistic_time_rate },
    { "custom",               custom_state }
};

//! Convert `IntegratedStateType` to `json`.
void to_json( json& jsonObject, const IntegratedStateType& integratedStateType );

//! Convert `json` to `IntegratedStateType`.
void from_json( const json& jsonObject, IntegratedStateType& integratedStateType );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_PROPAGATOR_H
