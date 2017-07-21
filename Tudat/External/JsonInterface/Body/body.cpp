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

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `BodySettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< BodySettings >& bodySettings )
{
    using namespace json_interface;
    using Keys = Keys::Body;

    // Initialise
    jsonObject = json( );

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
    if ( getValuePointer< json >( settings, keyTree + Keys::atmosphere ) )
    {
        bodySettings->atmosphereSettings = createAtmosphereSettings( settings, keyTree + Keys::atmosphere );
    }

    // Ephemeris
    if ( getValuePointer< json >( settings, keyTree + Keys::ephemeris ) )
    {
        bodySettings->ephemerisSettings = createEphemerisSettings( settings, keyTree + Keys::ephemeris );
    }

    // Gravity field
    if ( getValuePointer< json >( settings, keyTree + Keys::gravityField ) )
    {
        bodySettings->gravityFieldSettings = createGravityFieldSettings( settings, keyTree + Keys::gravityField );
    }

    // Rotation model
    if ( getValuePointer< json >( settings, keyTree + Keys::rotationModel ) )
    {
        bodySettings->rotationModelSettings = createRotationModelSettings( settings, keyTree + Keys::rotationModel );
    }

    // Shape model
    if ( getValuePointer< json >( settings, keyTree + Keys::shapeModel ) )
    {
        bodySettings->shapeModelSettings = createShapeModelSettings( settings, keyTree + Keys::shapeModel );
    }

    // Radiation pressure
    boost::shared_ptr< std::map< std::string, json > > jsonRadiationPressureInterfaces =
            getValuePointer< std::map< std::string, json > >( settings, keyTree + Keys::radiationPressure );
    if ( jsonRadiationPressureInterfaces )
    {
        std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;
        for ( auto entry : *jsonRadiationPressureInterfaces )
        {
            const std::string radiatingBody = entry.first;
            radiationPressureSettings[ radiatingBody ] = createRadiationPressureInterfaceSettings(
                        settings, radiatingBody, keyTree + Keys::radiationPressure, fallbackArea );
        }
        bodySettings->radiationPressureSettings = radiationPressureSettings;
    }

    // Aerodynamics
    if ( getValuePointer< json >( settings, keyTree + Keys::aerodynamics ) )
    {
        bodySettings->aerodynamicCoefficientSettings = createAerodynamicCoefficientSettings(
                    settings, keyTree + Keys::aerodynamics, fallbackArea );
    }
}

} // namespace json_interface

} // namespace tudat
