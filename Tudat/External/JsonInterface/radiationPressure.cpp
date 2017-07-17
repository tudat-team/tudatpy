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

#include "radiationPressure.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings )
{
    using namespace json_interface;

    // Initialise
    jsonObject = json( );

    // Get type
    jsonObject[ "type" ] = stringFromEnum( radiationPressureInterfaceSettings->getRadiationPressureType( ),
                                           radiationPressureTypes );

    /// CannonBallRadiationPressureInterfaceSettings
    boost::shared_ptr< CannonBallRadiationPressureInterfaceSettings > cannonBallRadiationPressureInterfaceSettings =
            boost::dynamic_pointer_cast< CannonBallRadiationPressureInterfaceSettings >(
                radiationPressureInterfaceSettings );
    if ( cannonBallRadiationPressureInterfaceSettings )
    {
        // Reference area
        jsonObject[ "referenceArea" ] = cannonBallRadiationPressureInterfaceSettings->getArea( );

        // Radiation pressure coefficient
        jsonObject[ "radiationPressureCoefficient" ] =
                cannonBallRadiationPressureInterfaceSettings->getRadiationPressureCoefficient( );

        // Occulting bodies
        jsonObject[ "ocultingBodies" ] =
                cannonBallRadiationPressureInterfaceSettings->getOccultingBodies( );
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `RadiationPressureInterfaceSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::RadiationPressureInterfaceSettings > createRadiationPressureInterfaceSettings(
        const json &settings, const std::string& sourceBodyName, const double& fallbackArea )
{
    using namespace simulation_setup;

    // Get radiation pressure coefficient type (cannonBall by default)
    const RadiationPressureType radiationPressureType = enumFromString(
                getValue< std::string >( settings, "type", "cannonBall" ), radiationPressureTypes );

    if ( radiationPressureType == cannon_ball )
    {
        // Get reference area (use fallback value if not NaN when referenceArea is not provided)
        const double referenceArea = isnan( fallbackArea ) ? getValue< double >( settings, "referenceArea" )
                                                           : getValue( settings, "referenceArea", fallbackArea );

        // Get radiation pressure coefficient
        const double radiationPressureCoefficients = getValue< double >( settings, "radiationPressureCoefficient" );

        // Get list of occulting bodies (empty by default, i.e. ignore eclipses)
        const std::vector< std::string > occultingBodies = getValue< std::vector< std::string > >(
                    settings, "ocultingBodies", { } );

        // Create and return settings
        return boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    sourceBodyName, referenceArea, radiationPressureCoefficients, occultingBodies );
    }
    else
    {
        throw std::runtime_error( stringFromEnum( radiationPressureType, radiationPressureTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interfaces

} // namespace tudat
