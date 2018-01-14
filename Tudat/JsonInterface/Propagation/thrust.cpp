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
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings )
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
        boost::shared_ptr< ThrustDirectionFromStateGuidanceSettings > directionFromStateGuidanceSettings
                = boost::dynamic_pointer_cast< ThrustDirectionFromStateGuidanceSettings >( directionSettings );
        assertNonNullPointer( directionFromStateGuidanceSettings );
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
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings )
{
    using namespace json_interface;
    using K = Keys::Propagator::Acceleration::Thrust::Direction;

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
        handleUnimplementedEnumValue( directionType, thrustDirectionTypes, unsupportedThrustDirectionTypes );
    }
}


// ThrustEngineSettings

//! Create a `json` object from a shared pointer to a `ThrustEngineSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ThrustEngineSettings >& magnitudeSettings )
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
        boost::shared_ptr< ConstantThrustEngineSettings > contantMagnitudeSettings
                = boost::dynamic_pointer_cast< ConstantThrustEngineSettings >( magnitudeSettings );
        assertNonNullPointer( contantMagnitudeSettings );
        jsonObject[ K::constantMagnitude ] = contantMagnitudeSettings->thrustMagnitude_;
        jsonObject[ K::specificImpulse ] = contantMagnitudeSettings->specificImpulse_;
        jsonObject[ K::bodyFixedDirection ] = contantMagnitudeSettings->bodyFixedThrustDirection_;
        return;
    }
    case from_engine_properties_thrust_magnitude:
    {
        boost::shared_ptr< FromBodyThrustEngineSettings > fromBodyMagnitudeSettings
                = boost::dynamic_pointer_cast< FromBodyThrustEngineSettings >( magnitudeSettings );
        assertNonNullPointer( fromBodyMagnitudeSettings );
        jsonObject[ K::useAllEngines ] = fromBodyMagnitudeSettings->useAllEngines_;
        return;
    }
    default:
        handleUnimplementedEnumValue( magnitudeType, thrustMagnitudeTypes, unsupportedThrustMagnitudeTypes );
    }
}

//! Create a shared pointer to a `ThrustEngineSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< ThrustEngineSettings >& magnitudeSettings )
{
    using namespace json_interface;
    using K = Keys::Propagator::Acceleration::Thrust::Magnitude;

    const ThrustMagnitudeTypes magnitudeType = getValue< ThrustMagnitudeTypes >( jsonObject, K::type );

    switch ( magnitudeType ) {
    case constant_thrust_magnitude:
    {
        ConstantThrustEngineSettings defaults( TUDAT_NAN, TUDAT_NAN );
        magnitudeSettings = boost::make_shared< ConstantThrustEngineSettings >(
                    getValue< double >( jsonObject, K::constantMagnitude ),
                    getValue< double >( jsonObject, K::specificImpulse ),
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
        handleUnimplementedEnumValue( magnitudeType, thrustMagnitudeTypes, unsupportedThrustMagnitudeTypes );
    }
}


// Thrust

//! Create a `json` object from a shared pointer to a `ThrustAccelerationSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ThrustAccelerationSettings >& thrustAccelerationSettings )
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
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< ThrustAccelerationSettings >& thrustAccelerationSettings )
{
    using namespace json_interface;
    using namespace interpolators;
    using K = Keys::Propagator::Acceleration::Thrust;

    if ( isDefined( jsonObject, K::dataInterpolation ) )
    {
        thrustAccelerationSettings = boost::make_shared< ThrustAccelerationSettings >(
                    getValue< boost::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > >(
                        jsonObject, K::dataInterpolation ),
                    getValue< double >( jsonObject, K::specificImpulse ),
                    getValue< ThrustFrames >( jsonObject, K::frame ),
                    getValue< std::string >( jsonObject, K::centralBody ) );
    }
    else
    {
        thrustAccelerationSettings = boost::make_shared< ThrustAccelerationSettings >(
                    getValue< boost::shared_ptr< ThrustDirectionGuidanceSettings > >( jsonObject, K::direction ),
                    getValue< boost::shared_ptr< ThrustEngineSettings > >( jsonObject, K::magnitude ) );
    }
}

} // namespace simulation_setup

} // namespace tudat
