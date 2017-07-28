/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/SimulationSetup/EnvironmentSetup/createGravityFieldVariations.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace gravitation
{

//! Map of `BodyDeformationTypes` supported by `json_interface`.
static std::map< std::string, BodyDeformationTypes > bodyDeformationTypes =
{
    { "basicSolidBody",     basic_solid_body },
    { "tabulatedVariation", tabulated_variation }
};

//! Convert `BodyDeformationTypes` to `json`.
void to_json( json& jsonObject, const BodyDeformationTypes& bodyDeformationType );

//! Convert `json` to `BodyDeformationTypes`.
void from_json( const json& jsonObject, BodyDeformationTypes& bodyDeformationType );

} // namespace gravitation

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GravityFieldVariationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< GravityFieldVariationSettings >& variationSettings );

//! Create a shared pointer to a `GravityFieldVariationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< GravityFieldVariationSettings >& variationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_GRAVITYFIELDVARIATION_H
