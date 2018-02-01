/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_REFERENCEFRAMES_H
#define TUDAT_JSONINTERFACE_REFERENCEFRAMES_H

#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace reference_frames
{

// AerodynamicsReferenceFrames

//! Map of `AerodynamicsReferenceFrames` string representations.
static std::map< AerodynamicsReferenceFrames, std::string > aerodynamicsReferenceFrames =
{
    { inertial_frame, "intertial" },
    { corotating_frame, "corotating" },
    { vertical_frame, "vertical" },
    { trajectory_frame, "trajectory" },
    { aerodynamic_frame, "aerodynamic" },
    { body_frame, "body" }
};

//! `AerodynamicsReferenceFrames` not supported by `json_interface`.
static std::vector< AerodynamicsReferenceFrames > unsupportedAerodynamicsReferenceFrames = { };

//! Convert `AerodynamicsReferenceFrames` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AerodynamicsReferenceFrames& aerodynamicsReferenceFrame )
{
    jsonObject = json_interface::stringFromEnum( aerodynamicsReferenceFrame, aerodynamicsReferenceFrames );
}

//! Convert `json` to `AerodynamicsReferenceFrames`.
inline void from_json( const nlohmann::json& jsonObject, AerodynamicsReferenceFrames& aerodynamicsReferenceFrame )
{
    aerodynamicsReferenceFrame =
            json_interface::enumFromString( jsonObject, aerodynamicsReferenceFrames );
}


// AerodynamicsReferenceFrameAngles

//! Map of `AerodynamicsReferenceFrameAngles` string representations.
static std::map< AerodynamicsReferenceFrameAngles, std::string > aerodynamicsReferenceFrameAngles =
{
    { latitude_angle, "latitude" },
    { longitude_angle, "longitude" },
    { heading_angle, "heading" },
    { flight_path_angle, "flightPath" },
    { angle_of_attack, "angleOfAttach" },
    { angle_of_sideslip, "sideslip" },
    { bank_angle, "bank" }
};

//! `AerodynamicsReferenceFrameAngles` not supported by `json_interface`.
static std::vector< AerodynamicsReferenceFrameAngles > unsupportedAerodynamicsReferenceFrameAngles = { };

//! Convert `AerodynamicsReferenceFrameAngles` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AerodynamicsReferenceFrameAngles& aerodynamicsReferenceFrameAngle )
{
    jsonObject = json_interface::stringFromEnum( aerodynamicsReferenceFrameAngle, aerodynamicsReferenceFrameAngles );
}

//! Convert `json` to `AerodynamicsReferenceFrameAngles`.
inline void from_json( const nlohmann::json& jsonObject, AerodynamicsReferenceFrameAngles& aerodynamicsReferenceFrameAngle )
{
    aerodynamicsReferenceFrameAngle =
            json_interface::enumFromString( jsonObject, aerodynamicsReferenceFrameAngles );
}

} // namespace reference_frames

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_REFERENCEFRAMES_H
