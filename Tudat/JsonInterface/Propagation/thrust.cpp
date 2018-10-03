/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Propagation/acceleration.h"
#include "Tudat/JsonInterface/Propagation/thrust.h"

#include "Tudat/JsonInterface/Mathematics/interpolation.h"

namespace tudat
{

namespace simulation_setup
{

// ThrustDirectionGuidanceSettings

//! Create a `json` object from a shared pointer to a `AccelerationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings )
{
    if ( ! directionSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Propagator::Acceleration::Thrust::Direction;

    const ThrustDirectionGuidanceTypes directionType = directionSettings->thrustDirectionType_;
    jsonObject[ K::type ] = directionType;
    jsonObject[ K::relativeBody ] = directionSettings->relativeBody_;

    switch ( directionType ) {
    case colinear_with_state_segment_thrust_direction:
    {
        std::shared_ptr< ThrustDirectionFromStateGuidanceSettings > directionFromStateGuidanceSettings
                = std::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( directionSettings );
        assertNonnullptrPointer( directionFromStateGuidanceSettings );
        jsonObject[ K::colinearWithVelocity ] = directionFromStateGuidanceSettings->isColinearWithVelocity_;
        jsonObject[ K::towardsRelativeBody ] = directionFromStateGuidanceSettings->directionIsOppositeToVector_;
        return;
    }
    case thrust_direction_from_existing_body_orientation:
    {
        return;
    }
    default:
        handleUnimplementedEnumValue( directionType, thrustDirectionTypes, unsupportedThrustDirectionTypes );
    }
}

//! Create a shared pointer to a `ThrustDirectionGuidanceSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings )
{
    using namespace json_interface;
    using K = Keys::Propagator::Acceleration::Thrust::Direction;

    const ThrustDirectionGuidanceTypes directionType = getValue< ThrustDirectionGuidanceTypes >( jsonObject, K::type );
    const std::string relativeBody = getValue< std::string >( jsonObject, K::relativeBody );

    switch ( directionType ) {
    case colinear_with_state_segment_thrust_direction:
    {
        directionSettings = std::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                    relativeBody,
                    getValue< bool >( jsonObject, K::colinearWithVelocity),
                    getValue< bool >( jsonObject, K::towardsRelativeBody ) );
        return;
    }
    case thrust_direction_from_existing_body_orientation:
    {
        directionSettings = std::make_shared< ThrustDirectionGuidanceSettings >( directionType, relativeBody );
        return;
    }
    default:
        handleUnimplementedEnumValue( directionType, thrustDirectionTypes, unsupportedThrustDirectionTypes );
    }
}


// ThrustMagnitudeSettings

//! Create a `json` object from a shared pointer to a `ThrustMagnitudeSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ThrustMagnitudeSettings >& magnitudeSettings )
{
    if ( ! magnitudeSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Propagator::Acceleration::Thrust::Magnitude;

    const ThrustMagnitudeTypes magnitudeType = magnitudeSettings->thrustMagnitudeGuidanceType_;
    jsonObject[ K::type ] = magnitudeType;
    assignIfNotEmpty( jsonObject, K::originID, magnitudeSettings->thrustOriginId_ );

    switch ( magnitudeType ) {
    case constant_thrust_magnitude:
    {
        std::shared_ptr< ConstantThrustMagnitudeSettings > contantMagnitudeSettings
                = std::dynamic_pointer_cast< ConstantThrustMagnitudeSettings >( magnitudeSettings );
        assertNonnullptrPointer( contantMagnitudeSettings );
        jsonObject[ K::constantMagnitude ] = contantMagnitudeSettings->thrustMagnitude_;
        jsonObject[ K::specificImpulse ] = contantMagnitudeSettings->specificImpulse_;
        jsonObject[ K::bodyFixedDirection ] = contantMagnitudeSettings->bodyFixedThrustDirection_;
        return;
    }
    case from_engine_properties_thrust_magnitude:
    {
        std::shared_ptr< FromBodyThrustMagnitudeSettings > fromBodyMagnitudeSettings
                = std::dynamic_pointer_cast< FromBodyThrustMagnitudeSettings >( magnitudeSettings );
        assertNonnullptrPointer( fromBodyMagnitudeSettings );
        jsonObject[ K::useAllEngines ] = fromBodyMagnitudeSettings->useAllEngines_;
        return;
    }
    default:
        handleUnimplementedEnumValue( magnitudeType, thrustMagnitudeTypes, unsupportedThrustMagnitudeTypes );
    }
}

//! Create a shared pointer to a `ThrustMagnitudeSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ThrustMagnitudeSettings >& magnitudeSettings )
{
    using namespace json_interface;
    using K = Keys::Propagator::Acceleration::Thrust::Magnitude;

    const ThrustMagnitudeTypes magnitudeType = getValue< ThrustMagnitudeTypes >( jsonObject, K::type );

    switch ( magnitudeType ) {
    case constant_thrust_magnitude:
    {
        ConstantThrustMagnitudeSettings defaults( TUDAT_NAN, TUDAT_NAN );
        magnitudeSettings = std::make_shared< ConstantThrustMagnitudeSettings >(
                    getValue< double >( jsonObject, K::constantMagnitude ),
                    getValue< double >( jsonObject, K::specificImpulse ),
                    getValue( jsonObject, K::bodyFixedDirection, defaults.bodyFixedThrustDirection_ ) );
        return;
    }
    case from_engine_properties_thrust_magnitude:
    {
        FromBodyThrustMagnitudeSettings defaults;
        magnitudeSettings = std::make_shared< FromBodyThrustMagnitudeSettings >(
                    getValue( jsonObject, K::useAllEngines, defaults.useAllEngines_ ),
                    getValue( jsonObject, K::originID, defaults.thrustOriginId_ ) );
        return;
    }
    default:
        handleUnimplementedEnumValue( magnitudeType, thrustMagnitudeTypes, unsupportedThrustMagnitudeTypes );
    }
}


// Thrust

//! Create a `json` object from a shared pointer to a `ThrustAccelerationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ThrustAccelerationSettings >& thrustAccelerationSettings )
{
    if ( ! thrustAccelerationSettings )
    {
        return;
    }
    using namespace json_interface;
    using namespace interpolators;
    using K = Keys::Propagator::Acceleration::Thrust;

    jsonObject[ Keys::Propagator::Acceleration::type ] = basic_astrodynamics::thrust_acceleration;

    if ( thrustAccelerationSettings->dataInterpolationSettings_ )
    {
        jsonObject[ K::dataInterpolation ] = thrustAccelerationSettings->dataInterpolationSettings_;
        jsonObject[ K::specificImpulse ] = thrustAccelerationSettings->constantSpecificImpulse_;
        jsonObject[ K::frame ] = thrustAccelerationSettings->thrustFrame_;
        assignIfNotEmpty( jsonObject, K::centralBody, thrustAccelerationSettings->centralBody_ );
    }
    else
    {
        jsonObject[ K::direction ] = thrustAccelerationSettings->thrustDirectionGuidanceSettings_;
        jsonObject[ K::magnitude ] = thrustAccelerationSettings->thrustMagnitudeSettings_;
    }
}

//! Create a shared pointer to a `ThrustAccelerationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ThrustAccelerationSettings >& thrustAccelerationSettings )
{
    using namespace json_interface;
    using namespace interpolators;
    using K = Keys::Propagator::Acceleration::Thrust;

    if ( isDefined( jsonObject, K::dataInterpolation ) )
    {
        thrustAccelerationSettings = std::make_shared< ThrustAccelerationSettings >(
                    getValue< std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > >(
                        jsonObject, K::dataInterpolation ),
                    getValue< double >( jsonObject, K::specificImpulse ),
                    getValue< ThrustFrames >( jsonObject, K::frame ),
                    getValue< std::string >( jsonObject, K::centralBody ) );
    }
    else
    {
        thrustAccelerationSettings = std::make_shared< ThrustAccelerationSettings >(
                    getValue< std::shared_ptr< ThrustDirectionGuidanceSettings > >( jsonObject, K::direction ),
                    getValue< std::shared_ptr< ThrustMagnitudeSettings > >( jsonObject, K::magnitude ) );
    }
}

} // namespace simulation_setup

} // namespace tudat
