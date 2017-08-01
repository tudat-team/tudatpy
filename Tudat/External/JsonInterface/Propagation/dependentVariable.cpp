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

#include "dependentVariable.h"

#include "acceleration.h"
#include "torque.h"
#include "referenceFrames.h"

namespace tudat
{

namespace propagators
{

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< SingleDependentVariableSaveSettings >& variableSettings )
{
    if ( ! variableSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Propagator::DependentVariable;

    const PropagationDependentVariables variableType = variableSettings->variableType_;
    jsonObject[ K::name ] = variableType;
    jsonObject[ K::associatedBody ] = variableSettings->associatedBody_;
    assignIfNotEmpty( jsonObject, K::secondaryBody, variableSettings->secondaryBody_ );

    switch ( variableType )
    {
    case single_acceleration_norm_dependent_variable:
    case single_acceleration_dependent_variable:
    {
        boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationVariableSettings =
                boost::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >( variableSettings );
        enforceNonNullPointer( accelerationVariableSettings );
        jsonObject[ K::accelerationType ] = accelerationVariableSettings->accelerationModeType_;
        return;
    }
    case single_torque_norm_dependent_variable:
    case single_torque_dependent_variable:
    {
        boost::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueVariableSettings =
                boost::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >( variableSettings );
        enforceNonNullPointer( torqueVariableSettings );
        jsonObject[ K::torqueType ] = torqueVariableSettings->torqueModeType_;
        return;
    }
    case intermediate_aerodynamic_rotation_matrix_variable:
    {
        boost::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings > aerodynamicRotationVariableSettings =
                boost::dynamic_pointer_cast< IntermediateAerodynamicRotationVariableSaveSettings >( variableSettings );
        enforceNonNullPointer( aerodynamicRotationVariableSettings );
        jsonObject[ K::baseFrame ] = aerodynamicRotationVariableSettings->baseFrame_;
        jsonObject[ K::targetFrame ] = aerodynamicRotationVariableSettings->targetFrame_;
        return;
    }
    case relative_body_aerodynamic_orientation_angle_variable:
    {
        boost::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > aerodynamicAngleVariableSettings =
                boost::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( variableSettings );
        enforceNonNullPointer( aerodynamicAngleVariableSettings );
        jsonObject[ K::angle ] = aerodynamicAngleVariableSettings->angle_;
        return;
    }
    default:
        return;
    }
}

//! Create a shared pointer to a `SingleDependentVariableSaveSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< SingleDependentVariableSaveSettings >& variableSettings )
{
    using namespace basic_astrodynamics;
    using namespace reference_frames;
    using namespace json_interface;
    using K = Keys::Propagator::DependentVariable;

    const PropagationDependentVariables variableName = getValue< PropagationDependentVariables >( jsonObject, K::name );
    const std::string associatedBody =
            getValue( jsonObject, K::bodyUndergoingTorque,
                      getValue( jsonObject, K::bodyUndergoingAcceleration,
                                getValue< std::string>( jsonObject, K::associatedBody ) ) );
    const std::string secondaryBody =
            getValue( jsonObject, K::bodyExertingTorque,
                      getValue( jsonObject, K::bodyExertingAcceleration,
                                getValue< std::string>( jsonObject, K::secondaryBody, "" ) ) );

    switch ( variableName )
    {
    case single_acceleration_norm_dependent_variable:
    case single_acceleration_dependent_variable:
    {
        variableSettings = boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    getValue< AvailableAcceleration >( jsonObject, K::accelerationType ),
                    associatedBody, secondaryBody, variableName == single_acceleration_norm_dependent_variable );
        return;
    }
    case single_torque_norm_dependent_variable:
    case single_torque_dependent_variable:
    {
        variableSettings = boost::make_shared< SingleTorqueDependentVariableSaveSettings >(
                    getValue< AvailableTorque >( jsonObject, K::torqueType ),
                    associatedBody, secondaryBody, variableName == single_torque_norm_dependent_variable );
        return;
    }
    case intermediate_aerodynamic_rotation_matrix_variable:
    {
        variableSettings = boost::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    associatedBody,
                    getValue< AerodynamicsReferenceFrames >( jsonObject, K::baseFrame ),
                    getValue< AerodynamicsReferenceFrames >( jsonObject, K::targetFrame ) );
        return;
    }
    case relative_body_aerodynamic_orientation_angle_variable:
    {
        variableSettings = boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    associatedBody, getValue< AerodynamicsReferenceFrameAngles >( jsonObject, K::angle ) );
        return;
    }
    default:
    {
        variableSettings = boost::make_shared< SingleDependentVariableSaveSettings >(
                    variableName, associatedBody, secondaryBody );
        return;
    }
    }
}

} // namespace propagators

} // namespace tudat
