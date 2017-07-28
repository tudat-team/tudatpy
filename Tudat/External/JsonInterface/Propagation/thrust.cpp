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

#include "thrust.h"

namespace tudat
{

namespace simulation_setup
{

/// DIRECTION

//! Convert `ThrustDirectionGuidanceTypes` to `json`.
void to_json( json& jsonObject, const ThrustDirectionGuidanceTypes& directionType )
{
    jsonObject = json_interface::stringFromEnum( directionType, thrustDirectionTypes );
}

//! Convert `json` to `ThrustDirectionGuidanceTypes`.
void from_json( const json& jsonObject, ThrustDirectionGuidanceTypes& directionType )
{
    directionType = json_interface::enumFromString( jsonObject.get< std::string >( ), thrustDirectionTypes );
}

//! Create a `json` object from a shared pointer to a `AccelerationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings )
{
    if ( ! directionSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Acceleration::ThrustDirection;

    const ThrustDirectionGuidanceTypes directionType = directionSettings->thrustDirectionType_;
    jsonObject[ K::type ] = directionType;
    jsonObject[ K::relativeBody ] = directionSettings->relativeBody_;

    switch ( directionType ) {
    case colinear_with_state_segment_thrust_direction:
    {
        boost::shared_ptr< ThrustDirectionFromStateGuidanceSettings > directionFromStateGuidanceSettings
                = boost::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( directionSettings );
        enforceNonNullPointer( directionFromStateGuidanceSettings );
        jsonObject[ K::colinearWithVelocity ] = directionFromStateGuidanceSettings->isColinearWithVelocity_;
        jsonObject[ K::towardsRelativeBody ] = directionFromStateGuidanceSettings->directionIsOppositeToVector_;
        return;
    }
    case thrust_direction_from_existing_body_orientation:
    {
        return;
    }
    default:
        jsonObject = handleUnimplementedEnumValueToJson( directionType, thrustDirectionTypes,
                                                         unsupportedThrustDirectionTypes );
    }
}

//! Create a shared pointer to a `ThrustDirectionGuidanceSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings )
{
    using namespace json_interface;
    using K = Keys::Acceleration::ThrustDirection;

    const ThrustDirectionGuidanceTypes directionType = getValue< ThrustDirectionGuidanceTypes >( jsonObject, K::type );
    const std::string relativeBody = getValue< std::string >( jsonObject, K::relativeBody );

    switch ( directionType ) {
    case colinear_with_state_segment_thrust_direction:
    {
        directionSettings = boost::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                    relativeBody,
                    getValue< bool >( jsonObject, K::colinearWithVelocity),
                    getValue< bool >( jsonObject, K::towardsRelativeBody ) );
        return;
    }
    case thrust_direction_from_existing_body_orientation:
    {
        directionSettings = boost::make_shared< ThrustDirectionGuidanceSettings >( directionType, relativeBody );
        return;
    }
    default:
        handleUnimplementedEnumValueFromJson( directionType, thrustDirectionTypes, unsupportedThrustDirectionTypes );
    }
}


/// MAGNITUDE

//! Convert `ThrustMagnitudeTypes` to `json`.
void to_json( json& jsonObject, const ThrustMagnitudeTypes& magnitudeType )
{
    jsonObject = json_interface::stringFromEnum( magnitudeType, thrustMagnitudeTypes );
}

//! Convert `json` to `ThrustMagnitudeTypes`.
void from_json( const json& jsonObject, ThrustMagnitudeTypes& magnitudeType )
{
    magnitudeType = json_interface::enumFromString( jsonObject.get< std::string >( ), thrustMagnitudeTypes );
}

//! Create a `json` object from a shared pointer to a `ThrustEngineSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ThrustEngineSettings >& magnitudeSettings )
{
    if ( ! magnitudeSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Acceleration::ThrustMagnitude;

    const ThrustMagnitudeTypes magnitudeType = magnitudeSettings->thrustMagnitudeGuidanceType_;
    jsonObject[ K::type ] = magnitudeType;
    jsonObject[ K::originID ] = magnitudeSettings->thrustOriginId_;

    switch ( magnitudeType ) {
    case constant_thrust_magnitude:
    {
        boost::shared_ptr< ConstantThrustEngineSettings > contantMagnitudeSettings
                = boost::dynamic_pointer_cast< ConstantThrustEngineSettings >( magnitudeSettings );
        enforceNonNullPointer( contantMagnitudeSettings );
        jsonObject[ K::constantMagnitude ] = contantMagnitudeSettings->thrustMagnitude_;
        jsonObject[ K::specificImpulse ] = contantMagnitudeSettings->specificImpulse_;
        jsonObject[ K::bodyFixedDirection ] = contantMagnitudeSettings->bodyFixedThrustDirection_;
        return;
    }
    case from_engine_properties_thrust_magnitude:
    {
        boost::shared_ptr< FromBodyThrustEngineSettings > fromBodyMagnitudeSettings
                = boost::dynamic_pointer_cast< FromBodyThrustEngineSettings >( magnitudeSettings );
        enforceNonNullPointer( fromBodyMagnitudeSettings );
        jsonObject[ K::useAllEngines ] = fromBodyMagnitudeSettings->useAllEngines_;
        return;
    }
    default:
        jsonObject = handleUnimplementedEnumValueToJson( magnitudeType, thrustMagnitudeTypes,
                                                         unsupportedThrustMagnitudeTypes );
    }
}

//! Create a shared pointer to a `ThrustEngineSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ThrustEngineSettings >& magnitudeSettings )
{
    using namespace json_interface;
    using K = Keys::Acceleration::ThrustMagnitude;

    const ThrustMagnitudeTypes magnitudeType = getValue< ThrustMagnitudeTypes >( jsonObject, K::type );

    switch ( magnitudeType ) {
    case constant_thrust_magnitude:
    {
        ConstantThrustEngineSettings defaults( TUDAT_NAN, TUDAT_NAN );
        magnitudeSettings = boost::make_shared< ConstantThrustEngineSettings >(
                    getNumeric< double >( jsonObject, K::constantMagnitude ),
                    getNumeric< double >( jsonObject, K::specificImpulse ),
                    getValue( jsonObject, K::bodyFixedDirection, defaults.bodyFixedThrustDirection_ ) );
        return;
    }
    case from_engine_properties_thrust_magnitude:
    {
        FromBodyThrustEngineSettings defaults;
        magnitudeSettings = boost::make_shared< FromBodyThrustEngineSettings >(
                    getValue( jsonObject, K::useAllEngines, defaults.useAllEngines_ ),
                    getValue( jsonObject, K::originID, defaults.thrustOriginId_ ) );
        return;
    }
    default:
        handleUnimplementedEnumValueFromJson( magnitudeType, thrustMagnitudeTypes, unsupportedThrustMagnitudeTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
