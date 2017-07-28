/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace ephemerides
{

//! Map of `ApproximatePlanetPositionsBase::BodiesWithEphemerisData` supported by `json_interface`.
static std::map< std::string, ApproximatePlanetPositionsBase::BodiesWithEphemerisData > bodiesWithEphemerisData =
{
    { "mercury",             ApproximatePlanetPositionsBase::mercury },
    { "venus",	             ApproximatePlanetPositionsBase::venus },
    { "earthMoonBarycenter", ApproximatePlanetPositionsBase::earthMoonBarycenter },
    { "mars",                ApproximatePlanetPositionsBase::mars },
    { "jupiter",             ApproximatePlanetPositionsBase::jupiter },
    { "saturn",              ApproximatePlanetPositionsBase::saturn },
    { "uranus",              ApproximatePlanetPositionsBase::uranus },
    { "neptune",             ApproximatePlanetPositionsBase::neptune },
    { "pluto",               ApproximatePlanetPositionsBase::pluto }
};

//! Convert `ApproximatePlanetPositionsBase::BodiesWithEphemerisData` to `json`.
void to_json( json& jsonObject, const ApproximatePlanetPositionsBase::BodiesWithEphemerisData& bodyWithEphemerisData );

//! Convert `json` to `ApproximatePlanetPositionsBase::BodiesWithEphemerisData`.
void from_json( const json& jsonObject, ApproximatePlanetPositionsBase::BodiesWithEphemerisData& bodyWithEphemerisData );

} // namespace ephemerides


namespace simulation_setup
{

//! Map of `EphemerisType`s supported by `json_interface`.
static std::map< std::string, EphemerisType > ephemerisTypes =
{
    { "approximatePlanetPositions", approximate_planet_positions },
    { "directSpice",                direct_spice_ephemeris },
    { "tabulated",                  tabulated_ephemeris },
    { "interpolatedSpice",		    interpolated_spice },
    { "constant",                   constant_ephemeris },
    { "kepler",                     kepler_ephemeris },
    // FIXME { "custom",                     custom_ephemeris }
};

//! Convert `EphemerisType` to `json`.
void to_json( json& jsonObject, const EphemerisType& ephemerisType );

//! Convert `json` to `EphemerisType`.
void from_json( const json& jsonObject, EphemerisType& ephemerisType );

//! Create a `json` object from a shared pointer to a `EphemerisSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< EphemerisSettings >& ephemerisSettings );

//! Create a shared pointer to a `EphemerisSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< EphemerisSettings >& ephemerisSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_EPHEMERIS_H
