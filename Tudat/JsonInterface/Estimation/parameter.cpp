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
#include "Tudat/JsonInterface/Estimation/parameter.h"
#include "Tudat/SimulationSetup/EstimationSetup/estimatableParameterSettings.h"

namespace tudat
{

namespace estimatable_parameters
{
// SingleDependentVariableSaveSettings

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( nlohmann::json& jsonObject,
              const boost::shared_ptr< EstimatableParameterSettings >& parameterSettings )
{
    if ( !parameterSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Parameter;

    const EstimatebleParametersEnum parameterType = parameterSettings->parameterType_.first;
    jsonObject[ K::parameterType ] = parameterType;
    jsonObject[ K::associatedBody ] = parameterSettings->parameterType_.second.first;

    switch ( parameterType )
    {
    case initial_body_state:
    {
        boost::shared_ptr< InitialTranslationalStateEstimatableParameterSettings< double > > stateEstimationSettings =
                boost::dynamic_pointer_cast< InitialTranslationalStateEstimatableParameterSettings< double > >(
                    parameterSettings );
        assertNonNullPointer( stateEstimationSettings );
        jsonObject[ K::initialStateValue ] = stateEstimationSettings->initialStateValue_;
        jsonObject[ K::centralBody ] = stateEstimationSettings->centralBody_;
        jsonObject[ K::frameOrientation ] = stateEstimationSettings->frameOrientation_;

        return;
    }
    case arc_wise_initial_body_state:
    {
        boost::shared_ptr< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > > stateEstimationSettings =
                boost::dynamic_pointer_cast<
                ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                    parameterSettings );
        assertNonNullPointer( stateEstimationSettings );
        jsonObject[ K::initialStateValue ] = stateEstimationSettings->initialStateValue_;
        jsonObject[ K::centralBody ] = stateEstimationSettings->centralBody_;
        jsonObject[ K::arcStartTimes ] = stateEstimationSettings->arcStartTimes_;
        jsonObject[ K::frameOrientation ] = stateEstimationSettings->frameOrientation_;

        return;
    }
    case gravitational_parameter:
    case radiation_pressure_coefficient:
    case constant_rotation_rate:
    case constant_drag_coefficient:
    case ppn_parameter_gamma:
    case ppn_parameter_beta:
    case equivalence_principle_lpi_violation_parameter:
    case direct_dissipation_tidal_time_lag:
    {
        boost::shared_ptr< DirectTidalTimeLagEstimatableParameterSettings > dissipationSettings =
                boost::dynamic_pointer_cast<
                DirectTidalTimeLagEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( dissipationSettings );
        jsonObject[ K::deformingBodies ] = dissipationSettings->deformingBodies_;
        return;
    }
    case constant_additive_observation_bias:
    {
        boost::shared_ptr< ConstantObservationBiasEstimatableParameterSettings > biasSettings =
                boost::dynamic_pointer_cast<
                ConstantObservationBiasEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( biasSettings );
        jsonObject[ K::observableType ] = biasSettings->observableType_;
        jsonObject[ K::linkEnds ] = biasSettings->linkEnds_;

    }
    case constant_relative_observation_bias:
    {
        boost::shared_ptr< ConstantObservationBiasEstimatableParameterSettings > biasSettings =
                boost::dynamic_pointer_cast<
                ConstantObservationBiasEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( biasSettings );
        jsonObject[ K::observableType ] = biasSettings->observableType_;
        jsonObject[ K::linkEnds ] = biasSettings->linkEnds_;
    }
    case arcwise_constant_additive_observation_bias:
    {
        boost::shared_ptr< ArcWiseConstantObservationBiasEstimatableParameterSettings > biasSettings =
                boost::dynamic_pointer_cast<
                ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( biasSettings );
        jsonObject[ K::observableType ] = biasSettings->observableType_;
        jsonObject[ K::linkEnds ] = biasSettings->linkEnds_;
        jsonObject[ K::arcStartTimes ] = biasSettings->arcStartTimes_;
        jsonObject[ K::referenceLinkEnd ] = biasSettings->linkEndForTime_;

    }
    case arcwise_constant_relative_observation_bias:
    {
        throw std::runtime_error( "JSON interface for bias estimation not yet imnplemented" );
    }
    case rotation_pole_position:
    case spherical_harmonics_cosine_coefficient_block:
    {
        throw std::runtime_error( "JSON interface for SH coefficient estimation not yet imnplemented" );
    }
    case spherical_harmonics_sine_coefficient_block:
    {
        throw std::runtime_error( "JSON interface for SH coefficient estimation not yet imnplemented" );
    }
    case ground_station_position:
    case empirical_acceleration_coefficients:
    {
        boost::shared_ptr< EmpiricalAccelerationEstimatableParameterSettings > empiricalAccelerationSettings =
                boost::dynamic_pointer_cast<
                EmpiricalAccelerationEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( empiricalAccelerationSettings );
        jsonObject[ K::centralBody ] = empiricalAccelerationSettings->centralBody_;
        jsonObject[ K::componentsToEstimate ] = empiricalAccelerationSettings->componentsToEstimate_;

        return;
    }
    case arc_wise_radiation_pressure_coefficient:
    {        
        boost::shared_ptr< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings > coefficientSettings =
                boost::dynamic_pointer_cast<
                ArcWiseRadiationPressureCoefficientEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( coefficientSettings );
        jsonObject[ K::arcStartTimes ] = coefficientSettings->arcStartTimeList_;

        return;
    }
    case arc_wise_empirical_acceleration_coefficients:
    {
        boost::shared_ptr< ArcWiseEmpiricalAccelerationEstimatableParameterSettings > empiricalAccelerationSettings =
                boost::dynamic_pointer_cast<
                ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( empiricalAccelerationSettings );
        jsonObject[ K::centralBody ] = empiricalAccelerationSettings->centralBody_;
        jsonObject[ K::componentsToEstimate ] = empiricalAccelerationSettings->componentsToEstimate_;
        jsonObject[ K::arcStartTimes ] = empiricalAccelerationSettings->arcStartTimeList_;

        return;
    }
    case full_degree_tidal_love_number:
    {

        boost::shared_ptr< FullDegreeTidalLoveNumberEstimatableParameterSettings > loveNumberSettings =
                boost::dynamic_pointer_cast<
                FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( loveNumberSettings );
        jsonObject[ K::deformingBodies ] = loveNumberSettings->deformingBodies_;
        jsonObject[ K::degree ] = loveNumberSettings->degree_;
        jsonObject[ K::useComplexValue ] = loveNumberSettings->useComplexValue_;

        return;
    }
    case single_degree_variable_tidal_love_number:
    {
        boost::shared_ptr< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings > loveNumberSettings =
                boost::dynamic_pointer_cast<
                SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( loveNumberSettings );
        jsonObject[ K::deformingBodies ] = loveNumberSettings->deformingBodies_;
        jsonObject[ K::degree ] = loveNumberSettings->degree_;
        jsonObject[ K::orders ] = loveNumberSettings->orders_;
        jsonObject[ K::useComplexValue ] = loveNumberSettings->useComplexValue_;

        return;
    }
    default:
    {
        assignIfNotEmpty( jsonObject, K::secondaryIdentifier, parameterSettings->parameterType_.second.second );
        return;
    }
    }
}

//! Create a shared pointer to a `SingleDependentVariableSaveSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject,
                boost::shared_ptr< EstimatableParameterSettings >& parameterSettings )
{
    using namespace basic_astrodynamics;
    using namespace reference_frames;
    using namespace json_interface;
    using K = Keys::Parameter;

    const EstimatebleParametersEnum parameterType =
            getValue< EstimatebleParametersEnum >( jsonObject, K::parameterType );
    const std::string bodyName = getValue< std::string>( jsonObject, K::associatedBody );

    switch ( parameterType )
    {
    case initial_body_state:
    {
        parameterSettings =
                boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                    bodyName,
                    getValue< Eigen::Matrix< double, 6, 1 > >( jsonObject, K::initialStateValue ),
                    getValue< std::string >( jsonObject, K::centralBody ),
                    getValue< std::string >( jsonObject, K::frameOrientation ) );

        return;
    }
    case arc_wise_initial_body_state:
    {
        parameterSettings =
                boost::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                    bodyName,
                    getValue< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( jsonObject, K::initialStateValue ),
                    getValue< std::vector< double > >( jsonObject, K::arcStartTimes ),
                    getValue< std::string >( jsonObject, K::centralBody ),
                    getValue< std::string >( jsonObject, K::frameOrientation ) );
        return;
    }
    case gravitational_parameter:
    case radiation_pressure_coefficient:
    case constant_rotation_rate:
    case constant_drag_coefficient:
    case ppn_parameter_gamma:
    case ppn_parameter_beta:
    case equivalence_principle_lpi_violation_parameter:
    case direct_dissipation_tidal_time_lag:
    {
        parameterSettings =
                boost::make_shared< DirectTidalTimeLagEstimatableParameterSettings >(
                    bodyName,
                    getValue< std::vector< std::string > >( jsonObject, K::deformingBodies ) );
        return;
    }
    case constant_additive_observation_bias:
    {
        throw std::runtime_error( "JSON interface for bias estimation not yet imnplemented" );
    }
    case constant_relative_observation_bias:
    {
        throw std::runtime_error( "JSON interface for bias estimation not yet imnplemented" );
    }
    case arcwise_constant_additive_observation_bias:
    {
        throw std::runtime_error( "JSON interface for bias estimation not yet imnplemented" );
    }
    case arcwise_constant_relative_observation_bias:
    {
        throw std::runtime_error( "JSON interface for bias estimation not yet imnplemented" );
    }
    case rotation_pole_position:
    case spherical_harmonics_cosine_coefficient_block:
    {
        throw std::runtime_error( "JSON interface for SH coefficient estimation not yet imnplemented" );
    }
    case spherical_harmonics_sine_coefficient_block:
    {
        throw std::runtime_error( "JSON interface for SH coefficient estimation not yet imnplemented" );
    }
    case ground_station_position:
    case empirical_acceleration_coefficients:
    {
        parameterSettings =
                boost::make_shared< EmpiricalAccelerationEstimatableParameterSettings >(
                    bodyName,
                    getValue< std::string >( jsonObject, K::centralBody ),
                    getValue< std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
                    std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > >(
                        jsonObject, K::componentsToEstimate ) );

        return;
    }
    case arc_wise_radiation_pressure_coefficient:
    {

    }
    case arc_wise_empirical_acceleration_coefficients:
    {

    }
    case full_degree_tidal_love_number:
    {

    }
    case single_degree_variable_tidal_love_number:
    {

    }
    default:
    {
        parameterSettings = boost::make_shared< EstimatableParameterSettings >(
                    bodyName,
                    parameterType,
                    getValue< std::string>( jsonObject, K::secondaryIdentifier, "" ) );
        return;
    }
    }
}

} // namespace propagators

} // namespace tudat
