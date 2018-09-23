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

#include "Tudat/JsonInterface/Environment/gravityFieldVariation.h"

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
    if ( !dependentVariableSettings )
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
        jsonObject[ K::accelerationType ] = accelerationVariableSettings->accelerationModelType_;
        jsonObject[ K::bodyExertingAcceleration ] = dependentVariableSettings->secondaryBody_;
        return;
    }
    case spherical_harmonic_acceleration_terms_dependent_variable:
    {
        boost::shared_ptr< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings > sphericalHarmonicsSettings =
                boost::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( sphericalHarmonicsSettings );
        jsonObject[ K::bodyExertingAcceleration ] = dependentVariableSettings->secondaryBody_;
        jsonObject[ K::componentIndices ] = sphericalHarmonicsSettings->componentIndices_;
        jsonObject[ K::componentIndex ] = dependentVariableSettings->componentIndex_;
        return;
    }
    case single_torque_norm_dependent_variable:
    case single_torque_dependent_variable:
    {
        boost::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueVariableSettings =
                boost::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( torqueVariableSettings );
        jsonObject[ K::torqueType ] = torqueVariableSettings->torqueModelType_;
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
    case single_gravity_field_variation_acceleration:
    {
        boost::shared_ptr< SingleVariationSphericalHarmonicAccelerationSaveSettings > variationalSphericalHarmonicsAccelerationSettings =
                boost::dynamic_pointer_cast< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( variationalSphericalHarmonicsAccelerationSettings );
        jsonObject[ K::bodyExertingAcceleration ] = variationalSphericalHarmonicsAccelerationSettings->secondaryBody_;
        jsonObject[ K::deformationType ] = variationalSphericalHarmonicsAccelerationSettings->deformationType_;
        jsonObject[ K::identifier ] = variationalSphericalHarmonicsAccelerationSettings->identifier_;
        return;
    }
    case single_gravity_field_variation_acceleration_terms:
    {
        boost::shared_ptr< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >
                variationalSphericalHarmonicsAccelerationTermsSettings =
                boost::dynamic_pointer_cast< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( variationalSphericalHarmonicsAccelerationTermsSettings );
        jsonObject[ K::bodyExertingAcceleration ] = variationalSphericalHarmonicsAccelerationTermsSettings->secondaryBody_;
        jsonObject[ K::componentIndices ] = variationalSphericalHarmonicsAccelerationTermsSettings->componentIndices_;
        jsonObject[ K::deformationType ] = variationalSphericalHarmonicsAccelerationTermsSettings->deformationType_;
        jsonObject[ K::identifier ] = variationalSphericalHarmonicsAccelerationTermsSettings->identifier_;
        return;
    }
    case acceleration_partial_wrt_body_translational_state:
    {
        boost::shared_ptr< AccelerationPartialWrtStateSaveSettings > accelerationPartialSettings =
                boost::dynamic_pointer_cast< AccelerationPartialWrtStateSaveSettings >(
                    dependentVariableSettings );
        assertNonNullPointer( accelerationPartialSettings );
        jsonObject[ K::bodyExertingAcceleration ] = accelerationPartialSettings->secondaryBody_;
        jsonObject[ K::componentIndices ] = accelerationPartialSettings->accelerationModelType_;
        jsonObject[ K::deformationType ] = accelerationPartialSettings->derivativeWrtBody_;
        jsonObject[ K::thirdBody ] = accelerationPartialSettings->thirdBody_;
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
    const std::string bodyName = getValue< std::string >( jsonObject, K::body );
    const int componentIndex = getValue< int >( jsonObject, K::componentIndex, -1 );

    switch ( dependentVariableType )
    {
    case single_acceleration_norm_dependent_variable:
    case single_acceleration_dependent_variable:
    {
        dependentVariableSettings = boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    getValue< AvailableAcceleration >( jsonObject, K::accelerationType ),
                    bodyName,
                    getValue< std::string >( jsonObject, K::bodyExertingAcceleration ),
                    dependentVariableType == single_acceleration_norm_dependent_variable,
                    componentIndex );
        return;
    }
    case spherical_harmonic_acceleration_terms_dependent_variable:
    {
        dependentVariableSettings = boost::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                    getValue< std::string >( jsonObject, K::body ),
                    getValue< std::string >( jsonObject, K::bodyExertingAcceleration ),
                    getValue< std::vector< std::pair< int, int > > >( jsonObject, K::componentIndices ),
                    componentIndex );
        return;
    }
    case single_torque_norm_dependent_variable:
    case single_torque_dependent_variable:
    {
        dependentVariableSettings = boost::make_shared< SingleTorqueDependentVariableSaveSettings >(
                    getValue< AvailableTorque >( jsonObject, K::torqueType ),
                    bodyName,
                    getValue< std::string >( jsonObject, K::bodyExertingTorque ),
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
                    getValue< AerodynamicsReferenceFrameAngles >( jsonObject, K::angle ),
                    getValue< std::string >( jsonObject, K::bodyExertingAcceleration, "" ) );
        return;
    }
    case single_gravity_field_variation_acceleration:
    {
        dependentVariableSettings = boost::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                    getValue< std::string >( jsonObject, K::body ),
                    getValue< std::string >( jsonObject, K::bodyExertingAcceleration ),
                    getValue< gravitation::BodyDeformationTypes >( jsonObject, K::deformationType ),
                    getValue< std::string >( jsonObject, K::identifier, "" ) );
        return;
    }
    case single_gravity_field_variation_acceleration_terms:
    {
        dependentVariableSettings = boost::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                    getValue< std::string >( jsonObject, K::body ),
                    getValue< std::string >( jsonObject, K::bodyExertingAcceleration ),
                    getValue< std::vector< std::pair< int, int > > >( jsonObject, K::componentIndices ),
                    getValue< gravitation::BodyDeformationTypes >( jsonObject, K::deformationType ),
                    getValue< std::string >( jsonObject, K::identifier, "" ) );
        return;
    }
    case acceleration_partial_wrt_body_translational_state:
    {
        dependentVariableSettings = boost::make_shared< AccelerationPartialWrtStateSaveSettings >(
                    getValue< std::string >( jsonObject, K::body ),
                    getValue< std::string >( jsonObject, K::bodyExertingAcceleration ),
                    getValue< basic_astrodynamics::AvailableAcceleration >( jsonObject, K::accelerationType ),
                    getValue< std::string >( jsonObject, K::derivativeWrtBody ),
                    getValue< std::string >( jsonObject, K::thirdBody, "" ) );
        return;
    }
    default:
    {
        dependentVariableSettings = boost::make_shared< SingleDependentVariableSaveSettings >(
                    dependentVariableType,
                    bodyName,
                    getValue< std::string >( jsonObject, K::relativeToBody, "" ),
                    componentIndex );
        return;
    }
    }
}

} // namespace propagators

} // namespace tudat
