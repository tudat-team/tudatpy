/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_EPHEMERIS_H
#define TUDAT_JSONINTERFACE_EPHEMERIS_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace ephemerides
{

//! Map of `ApproximatePlanetPositionsBase::BodiesWithEphemerisData`s string representations.
static std::map< ApproximatePlanetPositionsBase::BodiesWithEphemerisData, std::string >
bodiesWithEphemerisData =
{
    { ApproximatePlanetPositionsBase::mercury, "mercury" },
    { ApproximatePlanetPositionsBase::venus, "venus" },
    { ApproximatePlanetPositionsBase::earthMoonBarycenter, "earthMoonBarycenter" },
    { ApproximatePlanetPositionsBase::mars, "mars" },
    { ApproximatePlanetPositionsBase::jupiter, "jupiter" },
    { ApproximatePlanetPositionsBase::saturn, "saturn" },
    { ApproximatePlanetPositionsBase::uranus, "uranus" },
    { ApproximatePlanetPositionsBase::neptune, "neptune" },
    { ApproximatePlanetPositionsBase::pluto, "pluto" }
};

//! `ApproximatePlanetPositionsBase::BodiesWithEphemerisData` not supported by `json_interface`.
static std::vector< ApproximatePlanetPositionsBase::BodiesWithEphemerisData > unsupportedBodiesWithEphemerisData = { };

//! Convert `ApproximatePlanetPositionsBase::BodiesWithEphemerisData` to `json`.
inline void to_json( nlohmann::json& jsonObject,
                     const ApproximatePlanetPositionsBase::BodiesWithEphemerisData& bodyWithEphemerisData )
{
    jsonObject = json_interface::stringFromEnum( bodyWithEphemerisData, bodiesWithEphemerisData );
}

//! Convert `json` to `ApproximatePlanetPositionsBase::BodiesWithEphemerisData`.
inline void from_json( const nlohmann::json& jsonObject,
                       ApproximatePlanetPositionsBase::BodiesWithEphemerisData& bodyWithEphemerisData )
{
    bodyWithEphemerisData = json_interface::enumFromString( jsonObject, bodiesWithEphemerisData );
}

} // namespace ephemerides


namespace simulation_setup
{

//! Map of `EphemerisType`s string representations.
static std::map< EphemerisType, std::string > ephemerisTypes =
{
    { approximate_planet_positions, "approximatePlanetPositions" },
    { direct_spice_ephemeris, "directSpice" },
    { tabulated_ephemeris, "tabulated" },
    { interpolated_spice, "interpolatedSpice" },
    { constant_ephemeris, "constant" },
    { kepler_ephemeris, "kepler" },
    { custom_ephemeris, "custom" }
};

//! `EphemerisType` not supported by `json_interface`.
static std::vector< EphemerisType > unsupportedEphemerisTypes =
{
    custom_ephemeris
};

//! Convert `EphemerisType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const EphemerisType& ephemerisType )
{
    jsonObject = json_interface::stringFromEnum( ephemerisType, ephemerisTypes );
}

//! Convert `json` to `EphemerisType`.
inline void from_json( const nlohmann::json& jsonObject, EphemerisType& ephemerisType )
{
    ephemerisType = json_interface::enumFromString( jsonObject, ephemerisTypes );
}

//! Create a `json` object from a shared pointer to a `EphemerisSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< EphemerisSettings >& ephemerisSettings );

//! Create a shared pointer to a `EphemerisSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< EphemerisSettings >& ephemerisSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_EPHEMERIS_H
