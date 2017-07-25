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
    if ( bodySettings )
    {
        using namespace json_interface;
        using Keys = Keys::Body;

        // Constant mass
        const double constantMass = bodySettings->constantMass;
        if ( ! isnan( constantMass ) )
        {
            jsonObject[ Keys::mass ] = constantMass;
        }

        // Atmosphere
        if ( bodySettings->atmosphereSettings )
        {
            jsonObject[ Keys::atmosphere ] = bodySettings->atmosphereSettings;
        }

        // Ephemeris
        if ( bodySettings->ephemerisSettings )
        {
            jsonObject[ Keys::ephemeris ] = bodySettings->ephemerisSettings;
        }

        // Gravity field
        if ( bodySettings->gravityFieldSettings )
        {
            jsonObject[ Keys::gravityField ] = bodySettings->gravityFieldSettings;
        }

        // Rotation model
        if ( bodySettings->rotationModelSettings )
        {
            jsonObject[ Keys::rotationModel ] = bodySettings->rotationModelSettings;
        }

        // Shape model
        if ( bodySettings->shapeModelSettings )
        {
            jsonObject[ Keys::shapeModel ] = bodySettings->shapeModelSettings;
        }

        // Radiation pressure
        if ( ! bodySettings->radiationPressureSettings.empty( ) )
        {
            jsonObject[ Keys::radiationPressure ] = bodySettings->radiationPressureSettings;
        }

        // Aerodynamics
        if ( bodySettings->aerodynamicCoefficientSettings )
        {
            jsonObject[ Keys::aerodynamics ] = bodySettings->aerodynamicCoefficientSettings;
        }

        // Gravity field variations
        if ( ! bodySettings->gravityFieldVariationSettings.empty( ) )
        {
            jsonObject[ Keys::gravityFieldVariations ] = bodySettings->gravityFieldVariationSettings;
        }
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a `BodySettings` object with the settings from a `json` object.
boost::shared_ptr< simulation_setup::BodySettings > createBodySettings( const json& settings, const KeyTree& keyTree )
{
    using namespace simulation_setup;
    boost::shared_ptr< BodySettings > bodySettings = boost::make_shared< BodySettings >( );
    updateBodySettings( bodySettings, settings, keyTree );
    return bodySettings;
}

//! Update a `BodySettings` object with the settings from a `json` object.
void updateBodySettings( boost::shared_ptr< simulation_setup::BodySettings >& bodySettings,
                         const json& settings, const KeyTree& keyTree )
{
    using namespace simulation_setup;
    using Keys = Keys::Body;

    // Fallback reference area
    const double fallbackArea = getNumeric< double >( settings, keyTree + Keys::referenceArea, TUDAT_NAN, true );

    // Constant mass
    const boost::shared_ptr< double > mass = getNumericPointer< double >( settings, keyTree + Keys::mass );
    if ( mass )
    {
        bodySettings->constantMass = *mass ;
    }

    // Atmosphere
    if ( defined( settings, keyTree + Keys::atmosphere ) )
    {
        bodySettings->atmosphereSettings = createAtmosphereSettings( settings, keyTree + Keys::atmosphere );
    }

    // Ephemeris
    if ( defined( settings, keyTree + Keys::ephemeris ) )
    {
        bodySettings->ephemerisSettings = createEphemerisSettings( settings, keyTree + Keys::ephemeris );
    }

    // Gravity field
    if ( defined( settings, keyTree + Keys::gravityField ) )
    {
        bodySettings->gravityFieldSettings = createGravityFieldSettings( settings, keyTree + Keys::gravityField );
    }

    // Rotation model
    if ( defined( settings, keyTree + Keys::rotationModel ) )
    {
        bodySettings->rotationModelSettings = createRotationModelSettings( settings, keyTree + Keys::rotationModel );
    }

    // Shape model
    if ( defined( settings, keyTree + Keys::shapeModel ) )
    {
        bodySettings->shapeModelSettings = createShapeModelSettings( settings, keyTree + Keys::shapeModel );
    }

    // Radiation pressure
    if ( defined( settings, keyTree + Keys::radiationPressure ) )
    {
        std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;
        for ( auto entry : getValue< std::map< std::string, json > >( settings, keyTree + Keys::radiationPressure ) )
        {
            const std::string radiatingBody = entry.first;
            radiationPressureSettings[ radiatingBody ] = createRadiationPressureInterfaceSettings(
                        settings, radiatingBody, keyTree + Keys::radiationPressure, fallbackArea );
        }
        bodySettings->radiationPressureSettings = radiationPressureSettings;
    }

    // Aerodynamics
    if ( defined( settings, keyTree + Keys::aerodynamics ) )
    {
        bodySettings->aerodynamicCoefficientSettings = createAerodynamicCoefficientSettings(
                    settings, keyTree + Keys::aerodynamics, fallbackArea );
    }

    // Gravity field variations
    if ( defined( settings, keyTree + Keys::gravityFieldVariations ) )
    {
        std::vector< boost::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;
        for ( unsigned int i = 0;
              i < getValue< std::vector< json > >( settings, keyTree + Keys::gravityFieldVariations ).size( ); ++i )
        {
            gravityFieldVariationSettings.push_back(
                        createGravityFieldVariationSettings( settings, keyTree + Keys::gravityFieldVariations + i ) );
        }
        bodySettings->gravityFieldVariationSettings = gravityFieldVariationSettings;
    }
}

} // namespace json_interface

} // namespace tudat
