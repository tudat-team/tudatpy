/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "body.h"

#include "atmosphere.h"
#include "ephemeris.h"
#include "gravityField.h"
#include "rotationModel.h"
#include "shapeModel.h"
#include "radiationPressure.h"
#include "aerodynamics.h"
#include "gravityFieldVariation.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `BodySettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< BodySettings >& bodySettings )
{
    if ( ! bodySettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body;

    // Constant mass
    assignIfNotNaN( jsonObject, K::mass, bodySettings->constantMass );

    // Atmosphere
    if ( bodySettings->atmosphereSettings )
    {
        jsonObject[ K::atmosphere ] = bodySettings->atmosphereSettings;
    }

    // Ephemeris
    if ( bodySettings->ephemerisSettings )
    {
        jsonObject[ K::ephemeris ] = bodySettings->ephemerisSettings;
    }

    // Gravity field
    if ( bodySettings->gravityFieldSettings )
    {
        jsonObject[ K::gravityField ] = bodySettings->gravityFieldSettings;
    }

    // Rotation model
    if ( bodySettings->rotationModelSettings )
    {
        jsonObject[ K::rotationModel ] = bodySettings->rotationModelSettings;
    }

    // Shape model
    if ( bodySettings->shapeModelSettings )
    {
        jsonObject[ K::shapeModel ] = bodySettings->shapeModelSettings;
    }

    // Radiation pressure
    if ( ! bodySettings->radiationPressureSettings.empty( ) )
    {
        jsonObject[ K::radiationPressure ] = bodySettings->radiationPressureSettings;
    }

    // Aerodynamics
    if ( bodySettings->aerodynamicCoefficientSettings )
    {
        jsonObject[ K::aerodynamics ] = bodySettings->aerodynamicCoefficientSettings;
    }

    // Gravity field variations
    if ( ! bodySettings->gravityFieldVariationSettings.empty( ) )
    {
        jsonObject[ K::gravityFieldVariations ] = bodySettings->gravityFieldVariationSettings;
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a `BodySettings` object with the settings from a `json` object.
boost::shared_ptr< simulation_setup::BodySettings > createBodySettings( const json& jsonObject )
{
    using namespace simulation_setup;
    boost::shared_ptr< BodySettings > bodySettings = boost::make_shared< BodySettings >( );
    updateBodySettings( bodySettings, jsonObject );
    return bodySettings;
}

//! Update a `BodySettings` object with the settings from a `json` object.
void updateBodySettings( boost::shared_ptr< simulation_setup::BodySettings >& bodySettings, const json& jsonObject )
{
    using namespace simulation_setup;
    using K = Keys::Body;

    updateFromJSON( bodySettings->constantMass, jsonObject, K::mass, false );
    updateFromJSON( bodySettings->atmosphereSettings, jsonObject, K::atmosphere, false  );
    updateFromJSON( bodySettings->ephemerisSettings, jsonObject, K::ephemeris, false  );
    updateFromJSON( bodySettings->gravityFieldSettings, jsonObject, K::gravityField, false  );
    updateFromJSON( bodySettings->rotationModelSettings, jsonObject, K::rotationModel, false  );
    updateFromJSON( bodySettings->shapeModelSettings, jsonObject, K::shapeModel, false  );
    updateFromJSON( bodySettings->radiationPressureSettings, jsonObject, K::radiationPressure, false  );
    updateFromJSON( bodySettings->aerodynamicCoefficientSettings, jsonObject, K::aerodynamics, false  );
    updateFromJSON( bodySettings->gravityFieldVariationSettings, jsonObject, K::gravityFieldVariations, false  );

    /*
    // Gravity field variations
    if ( defined( jsonObject, K::gravityFieldVariations ) )
    {
        std::vector< boost::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;
        for ( unsigned int i = 0;
              i < getValue< std::vector< json > >( jsonObject, K::gravityFieldVariations ).size( ); ++i )
        {
            gravityFieldVariationSettings.push_back(
                        createGravityFieldVariationSettings( jsonObject, K::gravityFieldVariations / i ) );
        }
        bodySettings->gravityFieldVariationSettings = gravityFieldVariationSettings;
    }
    */
}

} // namespace json_interface

} // namespace tudat
