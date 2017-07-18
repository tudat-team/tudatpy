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
#include "aerodynamics.h"
#include "radiationPressure.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `BodySettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< BodySettings >& bodySettings )
{
    using namespace json_interface;

    // Initialise
    jsonObject = json( );

    // Atmosphere
    if ( bodySettings->atmosphereSettings )
    {
        jsonObject[ "atmosphere" ] = bodySettings->atmosphereSettings;
    }

    // Aerodynamics
    if ( bodySettings->aerodynamicCoefficientSettings )
    {
        jsonObject[ "aerodynamics" ] = bodySettings->aerodynamicCoefficientSettings;
    }

    // Radiation pressure
    if ( ! bodySettings->radiationPressureSettings.empty( ) )
    {
        jsonObject[ "radiationPressure" ] = bodySettings->radiationPressureSettings;
    }

    // Constant mass
    const double constantMass = bodySettings->constantMass;
    if ( constantMass == constantMass )
    {
        jsonObject[ "mass" ] = constantMass;
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Update a `BodySettings` object with the settings from a `json` object.
void updateBodySettings( std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > >& bodySettingsMap,
                         const std::string& bodyName, const json &settings )
{
    using namespace simulation_setup;

    // Create empty body settings if necessary
    if ( bodySettingsMap.count( bodyName ) == 0 )
    {
        bodySettingsMap[ bodyName ] = boost::make_shared< simulation_setup::BodySettings >( );
    }

    // Get reference to the body settings object to be updated
    boost::shared_ptr< simulation_setup::BodySettings > bodySettings = bodySettingsMap[ bodyName ];

    // Fallback reference area
    const double fallbackArea = getNumber< double >( settings, "referenceArea", TUDAT_NAN );

    /// Atmosphere
    boost::shared_ptr< json > jsonAtmosphereSettings = getValuePointer< json >( settings, "atmosphere" );
    if ( jsonAtmosphereSettings )
    {
        bodySettings->atmosphereSettings = createAtmosphereSettings( *jsonAtmosphereSettings );
    }

    /// Aerodynamics
    boost::shared_ptr< json > jsonAerodynamicCoefficientSettings = getValuePointer< json >( settings, "aerodynamics" );
    if ( jsonAerodynamicCoefficientSettings )
    {
        bodySettings->aerodynamicCoefficientSettings = createAerodynamicCoefficientSettings(
                    *jsonAerodynamicCoefficientSettings, fallbackArea );
    }

    /// Radiation pressure
    boost::shared_ptr< std::map< std::string, json > > jsonRadiationPressureInterfaces =
            getValuePointer< std::map< std::string, json > >( settings, "radiationPressure" );
    std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > > radiationPressureSettings;
    if ( jsonRadiationPressureInterfaces )
    {
        for ( auto entry : *jsonRadiationPressureInterfaces )
        {
            std::string radiatingBody = entry.first;
            json jsonRadiationPressureInterfaceSettings = entry.second;
            radiationPressureSettings[ radiatingBody ] = createRadiationPressureInterfaceSettings(
                        jsonRadiationPressureInterfaceSettings, radiatingBody, fallbackArea );
        }
        bodySettings->radiationPressureSettings = radiationPressureSettings;
    }

    /// Constant mass
    const boost::shared_ptr< double > mass = getNumberPointer< double >( settings, "mass" );
    if ( mass )
    {
        bodySettings->constantMass = *mass ;
    }
}

} // namespace json_interface

} // namespace tudat
