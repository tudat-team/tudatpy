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

#include "Tudat/JsonInterface/Estimation/observation.h"
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
    assignIfNotEmpty( jsonObject, K::secondaryIdentifier, parameterSettings->parameterType_.second.second );

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

        return;
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

        return;
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

        return;

    }
    case arcwise_constant_relative_observation_bias:
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

        return;
    }
    case spherical_harmonics_cosine_coefficient_block:
    {
        boost::shared_ptr< SphericalHarmonicEstimatableParameterSettings > coefficientSettings =
                boost::dynamic_pointer_cast<
                SphericalHarmonicEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( coefficientSettings );

        if( coefficientSettings->maximumDegree_ < 0 )
        {
            jsonObject[ K::coefficientIndices ] = coefficientSettings->blockIndices_;
        }
        else
        {
            jsonObject[ K::maximumDegree ] = coefficientSettings->maximumDegree_;
            jsonObject[ K::minimumDegree ] = coefficientSettings->minimumDegree_;
            jsonObject[ K::maximumOrder ] = coefficientSettings->maximumOrder_;
            jsonObject[ K::minimumOrder ] = coefficientSettings->minimumOrder_;

        }

        return;
    }
    case spherical_harmonics_sine_coefficient_block:
    {
        boost::shared_ptr< SphericalHarmonicEstimatableParameterSettings > coefficientSettings =
                boost::dynamic_pointer_cast<
                SphericalHarmonicEstimatableParameterSettings >(
                    parameterSettings );
        assertNonNullPointer( coefficientSettings );

        if( coefficientSettings->maximumDegree_ < 0 )
        {
            jsonObject[ K::coefficientIndices ] = coefficientSettings->blockIndices_;
        }
        else
        {
            jsonObject[ K::maximumDegree ] = coefficientSettings->maximumDegree_;
            jsonObject[ K::minimumDegree ] = coefficientSettings->minimumDegree_;
            jsonObject[ K::maximumOrder ] = coefficientSettings->maximumOrder_;
            jsonObject[ K::minimumOrder ] = coefficientSettings->minimumOrder_;

        }

        return;
    }
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
        parameterSettings = boost::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                    getValue< observation_models::LinkEnds >( jsonObject, K::linkEnds ),
                    getValue< observation_models::ObservableType >( jsonObject, K::observableType ), true );

        return;
    }
    case constant_relative_observation_bias:
    {
        parameterSettings = boost::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                    getValue< observation_models::LinkEnds >( jsonObject, K::linkEnds ),
                    getValue< observation_models::ObservableType >( jsonObject, K::observableType ), false );

        return;
    }
    case arcwise_constant_additive_observation_bias:
    {
        parameterSettings = boost::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                    getValue< observation_models::LinkEnds >( jsonObject, K::linkEnds ),
                    getValue< observation_models::ObservableType >( jsonObject, K::observableType ),
                    getValue< std::vector< double > >( jsonObject, K::arcStartTimes ),
                    getValue< observation_models::LinkEndType >( jsonObject, K::referenceLinkEnd ), true );

        return;
    }
    case arcwise_constant_relative_observation_bias:
    {
        parameterSettings = boost::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                    getValue< observation_models::LinkEnds >( jsonObject, K::linkEnds ),
                    getValue< observation_models::ObservableType >( jsonObject, K::observableType ),
                    getValue< std::vector< double > >( jsonObject, K::arcStartTimes ),
                    getValue< observation_models::LinkEndType >( jsonObject, K::referenceLinkEnd ), false );

        return;
    }
    case spherical_harmonics_cosine_coefficient_block:
    {
        bool parameterSet = false;

        try
        {
            parameterSettings =
                    boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                        getValue< int >( jsonObject, K::minimumDegree ),
                        getValue< int >( jsonObject, K::minimumOrder ),
                        getValue< int >( jsonObject, K::maximumDegree ),
                        getValue< int >( jsonObject, K::maximumOrder ),
                        bodyName,
                        spherical_harmonics_cosine_coefficient_block );

            parameterSet = true;
        }
        catch ( ... )
        {
            try
            {
                parameterSettings =
                        boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                            getValue< std::vector< std::pair< int, int > > >( jsonObject, K::coefficientIndices ),
                            bodyName,
                            spherical_harmonics_cosine_coefficient_block );
                if( parameterSet == true )
                {
                    std::cerr<<"Warning when reading spherical harmonic coefficient estimation settings from JSON, duplicate information detected"<<std::endl;
                }

                parameterSet = true;
            }
            catch ( ... ) { }
        }

        return;


    }
    case spherical_harmonics_sine_coefficient_block:
    {
        bool parameterSet = false;

        try
        {
            parameterSettings =
                    boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                        getValue< int >( jsonObject, K::minimumDegree ),
                        getValue< int >( jsonObject, K::minimumOrder ),
                        getValue< int >( jsonObject, K::maximumDegree ),
                        getValue< int >( jsonObject, K::maximumOrder ),
                        bodyName,
                        spherical_harmonics_sine_coefficient_block );

            parameterSet = true;
        }
        catch ( ... )
        {
            try
            {
                parameterSettings =
                        boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                            getValue< std::vector< std::pair< int, int > > >( jsonObject, K::coefficientIndices ),
                            bodyName,
                            spherical_harmonics_sine_coefficient_block );
                if( parameterSet == true )
                {
                    std::cerr<<"Warning when reading spherical harmonic coefficient estimation settings from JSON, duplicate information detected"<<std::endl;
                }

                parameterSet = true;
            }
            catch ( ... ) { }
        }

        return;
    }
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
        parameterSettings =
                boost::make_shared< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings >(
                    bodyName, getValue< std::vector< double > >( jsonObject, K::arcStartTimes ) );
        return;
    }
    case arc_wise_empirical_acceleration_coefficients:
    {

        parameterSettings =
                boost::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
                    bodyName,
                    getValue< std::string >( jsonObject, K::centralBody ),
                    getValue< std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
                    std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > >(
                        jsonObject, K::componentsToEstimate ),
                    getValue< std::vector< double > >( jsonObject, K::arcStartTimes ) );
        return;
    }
    case full_degree_tidal_love_number:
    {
        parameterSettings =
                boost::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >
                ( bodyName, getValue< int >( jsonObject, K::degree ),
                  getValue< std::vector< std::string > >( jsonObject, K::deformingBodies ),
                  getValue< bool >( jsonObject, K::useComplexValue ) );
        return;
    }
    case single_degree_variable_tidal_love_number:
    {
        parameterSettings =
                boost::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >
                ( bodyName, getValue< int >( jsonObject, K::degree ),
                  getValue< std::vector< int > >( jsonObject, K::orders ),
                  getValue< std::vector< std::string > >( jsonObject, K::deformingBodies ),
                  getValue< bool >( jsonObject, K::useComplexValue ) );
        return;
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
