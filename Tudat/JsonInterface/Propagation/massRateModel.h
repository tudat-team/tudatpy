/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_MASSRATEMODEL_H
#define TUDAT_JSONINTERFACE_MASSRATEMODEL_H

#include "Tudat/SimulationSetup/PropagationSetup/createMassRateModels.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Map of `AvailableMassRateModels` string representations.
static std::map< AvailableMassRateModels, std::string > massRateTypes =
{
    { undefined_mass_rate_model, "undefined" },
    { custom_mass_rate_model, "custom" },
    { from_thrust_mass_rate_model, "fromThrust" }
};

//! `AvailableMassRateModels` not supported by `json_interface`.
static std::vector< AvailableMassRateModels > unsupportedMassRateType =
{
    undefined_mass_rate_model,
    custom_mass_rate_model
};

//! Convert `AvailableMassRateModels` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AvailableMassRateModels& massRateType )
{
    jsonObject = json_interface::stringFromEnum( massRateType, massRateTypes );
}

//! Convert `json` to `AvailableMassRateModels`.
inline void from_json( const nlohmann::json& jsonObject, AvailableMassRateModels& massRateType )
{
    massRateType = json_interface::enumFromString( jsonObject, massRateTypes );
}

} // namespace basic_astrodynamics


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `MassRateModelSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< MassRateModelSettings >& massRateModelSettings );

//! Create a shared pointer to a `MassRateModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< MassRateModelSettings >& massRateModelSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_MASSRATEMODEL_H
