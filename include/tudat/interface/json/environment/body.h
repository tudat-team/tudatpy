/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/interface/json/environment/spice.h"

#include "tudat/interface/json/support/valueAccess.h"
#include "tudat/interface/json/support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `BodySettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< BodySettings >& bodySettings );

void to_json( nlohmann::json& jsonObject, const BodyListSettings& bodyListSettings );

} // namespace simulation_setup


namespace json_interface
{

//! Create a simulation_setup::BodySettings object with the settings from \p jsonObject.
/*!
 * @copybrief createBodySettings
 * \param jsonObject The `json` object containing the settings for one body.
 * \return Body settings object.
 */
std::shared_ptr< simulation_setup::BodySettings > createBodySettings( const nlohmann::json& jsonObject );

//! Update \p bodySettings with the settings from \p jsonObject.
/*!
 * @copybrief updateBodySettings
 * Does not change the values already defined in \p bodySettings that are not specified in \p jsonObject.
 * \param bodySettings Body settings object to be updated.
 * \param jsonObject The `json` object containing only the settings for one body.
 */
void updateBodySettings( std::shared_ptr< simulation_setup::BodySettings >& bodySettings, const nlohmann::json& jsonObject );

//! Update \p bodies and \p bodySettingsMap from \p jsonObject (using \p spiceSettings for default settings and
//! initial time from \p integratorSettings).
/*!
 * Update \p bodies and \p bodySettingsMap from \p jsonObject (using \p spiceSettings for default settings and
 * initial time from \p integratorSettings).
 * \param jsonObject The root `json` object containing all the relevant fields ("bodies" mandatory, "finalEpoch"
 * mandatory if Spice kernels should be preloaded, "globalFrameOrigin" and "globalFrameOrientation" optional).
 * \param bodies The named system of bodies to be updated (returned by reference).
 * \param bodySettingsMap The map of body settings created from \p jsonObject and used to create \p bodies
 * (returned by reference).
 * \param globalFrameOrigin Name of the global frame origin.
 * \param globalFrameOrientation Name of the global frame orientation.
 * \param spiceSettings The settings for Spice (nullptr if Spice not used).
 * \param integratorSettings The settings for the integrator (nullptr if Spice not used).
 * \throws std::runtime_error If any body is configured to be created using default settings and either
 * \p spiceSettings is `nullptr` or \p integratorSettings is `nullptr` and Spice is configured to preload kernels.
 */
template< typename TimeType = double >
void updateBodiesFromJSON(
        const nlohmann::json& jsonObject,
        simulation_setup::SystemOfBodies& bodies,
        simulation_setup::BodyListSettings& bodySettingsMap,
        const std::string globalFrameOrigin,
        const std::string globalFrameOrientation,
        const std::shared_ptr< SpiceSettings >& spiceSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& integratorSettings = nullptr )
{
    using namespace simulation_setup;

    bodySettingsMap.clear( );
    bodySettingsMap.resetFrames( globalFrameOrigin, globalFrameOrientation );

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
                                                              globalFrameOrigin,
                                                              globalFrameOrientation,
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
            std::shared_ptr< BodySettings > bodySettings = bodySettingsMap.at( bodyName );
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
            bodySettingsMap.addSettings( createBodySettings( jsonBodySettings ), bodyName );
        }
    }

    // Create bodies.
    bodies = createSystemOfBodies( bodySettingsMap );

    
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_BODY_H
