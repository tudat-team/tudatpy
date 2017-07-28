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

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace propagators
{

//! Map of `IntegratedStateType`s string representations.
static std::map< IntegratedStateType, std::string > integratedStateTypes =
{
    { hybrid, "hybrid" },
    { translational_state, "translational" },
    { rotational_state, "rotational" },
    { body_mass_state, "bodyMass" },
    { relativistic_time_rate, "relativisticTimeRate" },
    { custom_state, "custom" }
};

//! `IntegratedStateType`s not supported by `json_interface`.
static std::vector< IntegratedStateType > unsupportedIntegratedStateTypes = { };

//! Convert `IntegratedStateType` to `json`.
void to_json( json& jsonObject, const IntegratedStateType& integratedStateType );

//! Convert `json` to `IntegratedStateType`.
void from_json( const json& jsonObject, IntegratedStateType& integratedStateType );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_PROPAGATOR_H
