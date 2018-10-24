/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Environment/body.h"

#include "Tudat/JsonInterface/Environment/groundStations.h"
#include "Tudat/JsonInterface/Environment/atmosphere.h"
#include "Tudat/JsonInterface/Environment/ephemeris.h"
#include "Tudat/JsonInterface/Environment/gravityField.h"
#include "Tudat/JsonInterface/Environment/rotationModel.h"
#include "Tudat/JsonInterface/Environment/shapeModel.h"
#include "Tudat/JsonInterface/Environment/radiationPressure.h"
#include "Tudat/JsonInterface/Environment/aerodynamics.h"
#include "Tudat/JsonInterface/Environment/gravityFieldVariation.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `BodySettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< BodySettings >& bodySettings )
{
    if ( ! bodySettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body;

    assignIfNotNaN( jsonObject, K::mass, bodySettings->constantMass );
    assignIfNotnullptr( jsonObject, K::atmosphere, bodySettings->atmosphereSettings );
    assignIfNotnullptr( jsonObject, K::ephemeris, bodySettings->ephemerisSettings );
    assignIfNotnullptr( jsonObject, K::gravityField, bodySettings->gravityFieldSettings );
    assignIfNotnullptr( jsonObject, K::rotationModel, bodySettings->rotationModelSettings );
    assignIfNotnullptr( jsonObject, K::shapeModel, bodySettings->shapeModelSettings );
    assignIfNotEmpty( jsonObject, K::radiationPressure, bodySettings->radiationPressureSettings );
    assignIfNotnullptr( jsonObject, K::aerodynamics, bodySettings->aerodynamicCoefficientSettings );
    assignIfNotEmpty( jsonObject, K::gravityFieldVariation, bodySettings->gravityFieldVariationSettings );
    assignIfNotEmpty( jsonObject, K::groundStation, bodySettings->groundStationSettings );
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a simulation_setup::BodySettings object with the settings from \p jsonObject.
std::shared_ptr< simulation_setup::BodySettings > createBodySettings( const nlohmann::json& jsonObject )
{
    using namespace simulation_setup;
    std::shared_ptr< BodySettings > bodySettings = std::make_shared< BodySettings >( );
    updateBodySettings( bodySettings, jsonObject );
    return bodySettings;
}

//! Update \p bodySettings with the settings from \p jsonObject.
void updateBodySettings( std::shared_ptr< simulation_setup::BodySettings >& bodySettings, const nlohmann::json& jsonObject )
{
    using namespace simulation_setup;
    using K = Keys::Body;

    updateFromJSONIfDefined( bodySettings->constantMass, jsonObject, K::mass );
    updateFromJSONIfDefined( bodySettings->atmosphereSettings, jsonObject, K::atmosphere );
    updateFromJSONIfDefined( bodySettings->ephemerisSettings, jsonObject, K::ephemeris );
    updateFromJSONIfDefined( bodySettings->gravityFieldSettings, jsonObject, K::gravityField );
    updateFromJSONIfDefined( bodySettings->rotationModelSettings, jsonObject, K::rotationModel );
    updateFromJSONIfDefined( bodySettings->shapeModelSettings, jsonObject, K::shapeModel );
    updateFromJSONIfDefined( bodySettings->radiationPressureSettings, jsonObject, K::radiationPressure );
    updateFromJSONIfDefined( bodySettings->aerodynamicCoefficientSettings, jsonObject, K::aerodynamics );
    updateFromJSONIfDefined( bodySettings->gravityFieldVariationSettings, jsonObject, K::gravityFieldVariation );
    updateFromJSONIfDefined( bodySettings->groundStationSettings, jsonObject, K::groundStation );
}

} // namespace json_interface

} // namespace tudat
