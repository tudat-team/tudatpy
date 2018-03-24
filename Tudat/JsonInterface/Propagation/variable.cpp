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

#include "Tudat/JsonInterface/Propagation/variable.h"

#include "Tudat/JsonInterface/Propagation/acceleration.h"
#include "Tudat/JsonInterface/Propagation/torque.h"
#include "Tudat/JsonInterface/Propagation/referenceFrames.h"

namespace tudat
{

namespace propagators
{

// VariableSettings

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< VariableSettings >& variableSettings )
{
    if ( ! variableSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Variable;

    const VariableType variableType = variableSettings->variableType_;

    switch ( variableType )
    {
    case independentVariable:
    case cpuTimeVariable:
    case stateVariable:
    {
        jsonObject[ K::type ] = variableType;
        return;
    }
    case dependentVariable:
    {
        jsonObject = boost::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variableSettings );
        return;
    }
    default:
        handleUnimplementedEnumValue( variableType, variableTypes, unsupportedVariableTypes );
    }
}

//! Create a shared pointer to a `VariableSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< VariableSettings >& variableSettings )
{
    using namespace basic_astrodynamics;
    using namespace reference_frames;
    using namespace json_interface;
    using K = Keys::Variable;

    const VariableType variableType = getValue< VariableType >( jsonObject, K::type, dependentVariable );

    switch ( variableType )
    {
    case independentVariable:
    case cpuTimeVariable:
    case stateVariable:
    {
        variableSettings = boost::make_shared< VariableSettings >( variableType );
        return;
    }
    case dependentVariable:
    {
        variableSettings = getAs< boost::shared_ptr< SingleDependentVariableSaveSettings > >( jsonObject );
        return;
    }
    default:
        handleUnimplementedEnumValue( variableType, variableTypes, unsupportedVariableTypes );
    }
}


// SingleDependentVariableSaveSettings

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( nlohmann::json& jsonObject,
              const boost::shared_ptr< SingleDependentVariableSaveSettings >& dependentVariableSettings )
{
    if ( ! dependentVariableSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Variable;

    const PropagationDependentVariables dependentVariableType = dependentVariableSettings->dependentVariableType_;
    jsonObject[ K::dependentVariableType ] = dependentVariableType;
    jsonObject[ K::body ] = dependentVariableSettings->associatedBody_;
    if ( dependentVariableSettings->componentIndex_ >= 0 )
    {
        jsonObject[ K::componentIndex ] = dependentVariableSettings->componentIndex_;
    }

    switch ( dependentVariableType )
    {
    case single_acceleration_norm_dependent_variable:
    case single_acceleration_dependent_variable:
    {
        boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationVariableSettings =
                boost::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( accelerationVariableSettings );
        jsonObject[ K::accelerationType ] = accelerationVariableSettings->accelerationModeType_;
        jsonObject[ K::bodyExertingAcceleration ] = dependentVariableSettings->secondaryBody_;
        return;
    }
    case single_torque_norm_dependent_variable:
    case single_torque_dependent_variable:
    {
        boost::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueVariableSettings =
                boost::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( torqueVariableSettings );
        jsonObject[ K::torqueType ] = torqueVariableSettings->torqueModeType_;
        jsonObject[ K::bodyExertingTorque ] = dependentVariableSettings->secondaryBody_;
        return;
    }
    case intermediate_aerodynamic_rotation_matrix_variable:
    {
        boost::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings > aerodynamicRotationVariableSettings =
                boost::dynamic_pointer_cast< IntermediateAerodynamicRotationVariableSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( aerodynamicRotationVariableSettings );
        jsonObject[ K::baseFrame ] = aerodynamicRotationVariableSettings->baseFrame_;
        jsonObject[ K::targetFrame ] = aerodynamicRotationVariableSettings->targetFrame_;
        return;
    }
    case relative_body_aerodynamic_orientation_angle_variable:
    {
        boost::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > aerodynamicAngleVariableSettings =
                boost::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( aerodynamicAngleVariableSettings );
        jsonObject[ K::angle ] = aerodynamicAngleVariableSettings->angle_;
        return;
    }
    default:
    {
        assignIfNotEmpty( jsonObject, K::relativeToBody, dependentVariableSettings->secondaryBody_ );
        return;
    }
    }
}

//! Create a shared pointer to a `SingleDependentVariableSaveSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject,
                boost::shared_ptr< SingleDependentVariableSaveSettings >& dependentVariableSettings )
{
    using namespace basic_astrodynamics;
    using namespace reference_frames;
    using namespace json_interface;
    using K = Keys::Variable;

    const PropagationDependentVariables dependentVariableType =
            getValue< PropagationDependentVariables >( jsonObject, K::dependentVariableType );
    const std::string bodyName = getValue< std::string>( jsonObject, K::body );
    const int componentIndex = getValue< int >( jsonObject, K::componentIndex, -1 );

    switch ( dependentVariableType )
    {
    case single_acceleration_norm_dependent_variable:
    case single_acceleration_dependent_variable:
    {
        dependentVariableSettings = boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    getValue< AvailableAcceleration >( jsonObject, K::accelerationType ),
                    bodyName,
                    getValue< std::string>( jsonObject, K::bodyExertingAcceleration ),
                    dependentVariableType == single_acceleration_norm_dependent_variable,
                    componentIndex );
        return;
    }
    case single_torque_norm_dependent_variable:
    case single_torque_dependent_variable:
    {
        dependentVariableSettings = boost::make_shared< SingleTorqueDependentVariableSaveSettings >(
                    getValue< AvailableTorque >( jsonObject, K::torqueType ),
                    bodyName,
                    getValue< std::string>( jsonObject, K::bodyExertingTorque ),
                    dependentVariableType == single_torque_norm_dependent_variable,
                    componentIndex );
        return;
    }
    case intermediate_aerodynamic_rotation_matrix_variable:
    {
        dependentVariableSettings = boost::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    bodyName,
                    getValue< AerodynamicsReferenceFrames >( jsonObject, K::baseFrame ),
                    getValue< AerodynamicsReferenceFrames >( jsonObject, K::targetFrame ),
                    componentIndex );
        return;
    }
    case relative_body_aerodynamic_orientation_angle_variable:
    {
        dependentVariableSettings = boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    bodyName,
                    getValue< AerodynamicsReferenceFrameAngles >( jsonObject, K::angle ) );
        return;
    }
    default:
    {
        dependentVariableSettings = boost::make_shared< SingleDependentVariableSaveSettings >(
                    dependentVariableType,
                    bodyName,
                    getValue< std::string>( jsonObject, K::relativeToBody, "" ),
                    componentIndex );
        return;
    }
    }
}

} // namespace propagators

} // namespace tudat
