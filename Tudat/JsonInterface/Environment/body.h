/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_BODY_H
#define TUDAT_JSONINTERFACE_BODY_H

#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/JsonInterface/Environment/spice.h"

#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `BodySettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< BodySettings >& bodySettings );

} // namespace simulation_setup


namespace json_interface
{

//! Create a simulation_setup::BodySettings object with the settings from \p jsonObject.
/*!
 * @copybrief createBodySettings
 * \param jsonObject The `json` object containing the settings for one body.
 * \return Body settings object.
 */
boost::shared_ptr< simulation_setup::BodySettings > createBodySettings( const nlohmann::json& jsonObject );

//! Update \p bodySettings with the settings from \p jsonObject.
/*!
 * @copybrief updateBodySettings
 * Does not change the values already defined in \p bodySettings that are not specified in \p jsonObject.
 * \param bodySettings Body settings object to be updated.
 * \param jsonObject The `json` object containing only the settings for one body.
 */
void updateBodySettings( boost::shared_ptr< simulation_setup::BodySettings >& bodySettings, const nlohmann::json& jsonObject );

//! Update \p bodyMap and \p bodySettingsMap from \p jsonObject (using \p spiceSettings for default settings and
//! initial time from \p integratorSettings).
/*!
 * Update \p bodyMap and \p bodySettingsMap from \p jsonObject (using \p spiceSettings for default settings and
 * initial time from \p integratorSettings).
 * \param jsonObject The root `json` object containing all the relevant fields ("bodies" mandatory, "finalEpoch"
 * mandatory if Spice kernels should be preloaded, "globalFrameOrigin" and "globalFrameOrientation" optional).
 * \param bodyMap The named body map to be updated (returned by reference).
 * \param bodySettingsMap The map of body settings created from \p jsonObject and used to create \p bodyMap
 * (returned by reference).
 * \param globalFrameOrigin Name of the global frame origin.
 * \param globalFrameOrientation Name of the global frame orientation.
 * \param spiceSettings The settings for Spice (NULL if Spice not used).
 * \param integratorSettings The settings for the integrator (NULL if Spice not used).
 * \throws std::runtime_error If any body is configured to be created using default settings and either
 * \p spiceSettings is `NULL` or \p integratorSettings is `NULL` and Spice is configured to preload kernels.
 */
template< typename TimeType >
void updateBodiesFromJSON(
        const nlohmann::json& jsonObject,
        simulation_setup::NamedBodyMap& bodyMap,
        std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > >& bodySettingsMap,
        const std::string globalFrameOrigin,
        const std::string globalFrameOrientation,
        const boost::shared_ptr< SpiceSettings >& spiceSettings,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace simulation_setup;

    bodySettingsMap.clear( );

    std::map< std::string, nlohmann::json > jsonBodySettingsMap =
            getValue< std::map< std::string, nlohmann::json > >( jsonObject, Keys::bodies );

    std::vector< std::string > defaultBodyNames;
    for ( auto entry : jsonBodySettingsMap )
    {
        const std::string bodyName = entry.first;
        if ( getValue( jsonObject, Keys::bodies / bodyName / Keys::Body::useDefaultSettings, false ) )
        {
            defaultBodyNames.push_back( bodyName );
        }
    }

    // Create map with default body settings.
    if ( ! defaultBodyNames.empty( ) )
    {
        if ( spiceSettings )
        {
            if ( spiceSettings->preloadEphemeris_ )
            {
                if ( integratorSettings )
                {
                    const TimeType initialEpoch = integratorSettings->initialTime_;
                    const TimeType finalEpoch = getValue< TimeType >( jsonObject, Keys::finalEpoch );
                    const TimeType earliestInterpolationEpoch = std::min( initialEpoch, finalEpoch ) -
                            std::fabs( spiceSettings->getInitialOffset( ) );
                    const TimeType latestInterpolationEpoch = std::max( initialEpoch, finalEpoch ) +
                            std::fabs( spiceSettings->getFinalOffset( ) );
                    bodySettingsMap = getDefaultBodySettings( defaultBodyNames,
                                                              earliestInterpolationEpoch,
                                                              latestInterpolationEpoch,
                                                              spiceSettings->interpolationStep_ );
                }
                else
                {
                    throw std::runtime_error(
                                "Could not get default bodies settings because the initial time is not known. "
                                "Provide a valid IntegratorSettings object or turn Spice's preloadEphemeris off." );
                }
            }
            else
            {
                bodySettingsMap = getDefaultBodySettings( defaultBodyNames );
            }
        }
        else
        {
            bodySettingsMap = getDefaultBodySettings( defaultBodyNames );
        }
    }

    // Get body settings from JSON.
    for ( auto entry : jsonBodySettingsMap )
    {
        const std::string bodyName = entry.first;
        const nlohmann::json jsonBodySettings = jsonBodySettingsMap.at( bodyName );
        if ( bodySettingsMap.count( bodyName ) )
        {
            boost::shared_ptr< BodySettings >& bodySettings = bodySettingsMap[ bodyName ];
            // Reset ephemeris and rotational models frames.
            if ( bodySettings->ephemerisSettings )
            {
                bodySettings->ephemerisSettings->resetFrameOrientation( globalFrameOrientation );
            }
            if ( bodySettings->rotationModelSettings )
            {
                bodySettings->rotationModelSettings->resetOriginalFrame( globalFrameOrientation );
            }
            // Update body settings from JSON.
            updateBodySettings( bodySettings, jsonBodySettings );
        }
        else
        {
            // Create body settings from JSON.
            bodySettingsMap[ bodyName ] = createBodySettings( jsonBodySettings );
        }
    }

    // Create bodies.
    bodyMap = createBodies( bodySettingsMap );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, globalFrameOrigin, globalFrameOrientation );
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_BODY_H
