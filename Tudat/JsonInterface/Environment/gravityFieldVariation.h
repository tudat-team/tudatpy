/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_GRAVITYFIELDVARIATION_H
#define TUDAT_JSONINTERFACE_GRAVITYFIELDVARIATION_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityFieldVariations.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace gravitation
{

//! Map of `BodyDeformationTypes` string representations.
static std::map< BodyDeformationTypes, std::string > bodyDeformationTypes =
{
    { basic_solid_body, "basicSolidBody" },
    { tabulated_variation, "tabulatedVariation" }
};

//! `BodyDeformationTypes` not supported by `json_interface`.
static std::vector< BodyDeformationTypes > unsupportedBodyDeformationTypes =
{
};

//! Convert `BodyDeformationTypes` to `json`.
inline void to_json( nlohmann::json& jsonObject, const BodyDeformationTypes& bodyDeformationType )
{
    jsonObject = json_interface::stringFromEnum( bodyDeformationType, bodyDeformationTypes );
}

//! Convert `json` to `BodyDeformationTypes`.
inline void from_json( const nlohmann::json& jsonObject, BodyDeformationTypes& bodyDeformationType )
{
    bodyDeformationType = json_interface::enumFromString( jsonObject, bodyDeformationTypes );
}

} // namespace gravitation

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GravityFieldVariationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< GravityFieldVariationSettings >& variationSettings );

//! Create a shared pointer to a `GravityFieldVariationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< GravityFieldVariationSettings >& variationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_GRAVITYFIELDVARIATION_H
