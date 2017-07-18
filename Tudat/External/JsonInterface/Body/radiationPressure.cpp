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

//! Convert `RadiationPressureType` to `json`.
void to_json( json& jsonObject, const RadiationPressureType& radiationPressureType )
{
    jsonObject = json( json_interface::stringFromEnum( radiationPressureType, radiationPressureTypes ) );
}

//! Convert `json` to `RadiationPressureType`.
void from_json( const json& jsonObject, RadiationPressureType& radiationPressureType )
{
    radiationPressureType = json_interface::enumFromString( jsonObject.get< std::string >( ), radiationPressureTypes );
}

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings )
{
    using namespace json_interface;
    using Keys = Keys::Body::RadiationPressure;

    // Initialise
    jsonObject = json( );

    // Get type
    jsonObject[ Keys::type ] = radiationPressureInterfaceSettings->getRadiationPressureType( );

    /// CannonBallRadiationPressureInterfaceSettings
    boost::shared_ptr< CannonBallRadiationPressureInterfaceSettings > cannonBallRadiationPressureInterfaceSettings =
            boost::dynamic_pointer_cast< CannonBallRadiationPressureInterfaceSettings >(
                radiationPressureInterfaceSettings );
    if ( cannonBallRadiationPressureInterfaceSettings )
    {
        // Reference area
        jsonObject[ Keys::referenceArea ] = cannonBallRadiationPressureInterfaceSettings->getArea( );

        // Radiation pressure coefficient
        jsonObject[ Keys::radiationPressureCoefficient ] =
                cannonBallRadiationPressureInterfaceSettings->getRadiationPressureCoefficient( );

        // Occulting bodies
        jsonObject[ Keys::ocultingBodies ] =
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
    using Keys = Keys::Body::RadiationPressure;

    // Get radiation pressure coefficient type (cannonBall by default)
    const RadiationPressureType radiationPressureType = getValue( settings, Keys::type, cannon_ball );

    if ( radiationPressureType == cannon_ball )
    {
        // Get reference area (use fallback value if not NaN when referenceArea is not provided)
        const double referenceArea = isnan( fallbackArea ) ? getValue< double >( settings, Keys::referenceArea )
                                                           : getValue( settings, Keys::referenceArea, fallbackArea );

        // Get radiation pressure coefficient
        const double radiationPressureCoefficient = getValue< double >( settings, Keys::radiationPressureCoefficient );

        // Create settings
        CannonBallRadiationPressureInterfaceSettings cannonBallRadiationPressureInterfaceSettings(
                    sourceBodyName, referenceArea, radiationPressureCoefficient );

        // Get list of occulting bodies
        const auto occultingBodies = getValuePointer< std::vector< std::string > >( settings, Keys::ocultingBodies );
        if ( occultingBodies )
        {
            cannonBallRadiationPressureInterfaceSettings.occultingBodies_ = *occultingBodies;
        }

        // Return shared pointer
        return boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    cannonBallRadiationPressureInterfaceSettings );
    }
    else
    {
        throw std::runtime_error( stringFromEnum( radiationPressureType, radiationPressureTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
