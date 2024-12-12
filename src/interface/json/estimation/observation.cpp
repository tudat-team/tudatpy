/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/interface/json/estimation/observation.h"

#include "tudat/interface/json/support/utilities.h"

namespace tudat
{

namespace observation_models
{

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< ObservationSettings >& observationSettings )
{
    if( !observationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Observation;

    const ObservableType observableType = observationSettings->observableType_;
    jsonObject[ K::observableType ] = observableType;
    assignIfNotEmpty( jsonObject,
                      K::lightTimeCorrectionSettingsList,
                      observationSettings->lightTimeCorrectionsList_ );
    if( observationSettings->biasSettings_ != nullptr )
    {
        jsonObject[ K::biasSettings ] = observationSettings->biasSettings_;
    }

    switch( observableType )
    {
        case one_way_range:
        case angular_position:
        case relative_angular_position:
        case position_observable:
        case relative_position_observable:
        case one_way_doppler: {
            std::shared_ptr< OneWayDopplerObservationSettings > oneWayDopplerObservationSettings =
                    std::dynamic_pointer_cast< OneWayDopplerObservationSettings >(
                            observationSettings );
            if( !( oneWayDopplerObservationSettings == nullptr ) )
            {
                if( oneWayDopplerObservationSettings->transmitterProperTimeRateSettings_ !=
                    nullptr )
                {
                    jsonObject[ K::transmitterProperTimeRateSettings ] =
                            oneWayDopplerObservationSettings->transmitterProperTimeRateSettings_;
                }

                if( oneWayDopplerObservationSettings->receiverProperTimeRateSettings_ != nullptr )
                {
                    jsonObject[ K::receiverProperTimeRateSettings ] =
                            oneWayDopplerObservationSettings->receiverProperTimeRateSettings_;
                }
            }
            return;
        }
        case one_way_differenced_range: {
            std::shared_ptr< OneWayDifferencedRangeRateObservationSettings >
                    oneWayDifferencedRangeObservationSettings = std::dynamic_pointer_cast<
                            OneWayDifferencedRangeRateObservationSettings >( observationSettings );
            assertNonnullptrPointer( oneWayDifferencedRangeObservationSettings );

            std::cerr << "Warning, assuming constant integration time of 1-way differenced range "
                         "when parsing to JSON, correctness not ensured"
                      << std::endl;
            jsonObject[ K::constantIntegrationTime ] =
                    oneWayDifferencedRangeObservationSettings->integrationTimeFunction_( 0.0 );
            return;
        }
        case n_way_range: {
            std::shared_ptr< NWayRangeObservationSettings > nWayRangeObservationSettings =
                    std::dynamic_pointer_cast< NWayRangeObservationSettings >(
                            observationSettings );
            if( !( nWayRangeObservationSettings == nullptr ) )
            {
                jsonObject[ K::oneWayRangeObservationSettings ] =
                        nWayRangeObservationSettings->oneWayRangeObservationSettings_;

                if( !( nWayRangeObservationSettings->retransmissionTimesFunction_ ) )
                {
                    std::vector< double > retransmissionTimes =
                            nWayRangeObservationSettings->retransmissionTimesFunction_( 0.0 );
                    std::cerr << "Warning, assuming constant delay times of n-way range when "
                                 "parsing to JSON, correctness not ensured"
                              << std::endl;
                    jsonObject[ K::retransmissionTimes ] = retransmissionTimes;
                }
            }
            return;
        }
        case two_way_doppler: {
            std::shared_ptr< TwoWayDopplerObservationSettings > twoWayDopplerObservationSettings =
                    std::dynamic_pointer_cast< TwoWayDopplerObservationSettings >(
                            observationSettings );
            if( !( twoWayDopplerObservationSettings == nullptr ) )
            {
                jsonObject[ K::uplinkOneWayDopplerSettings ] =
                        twoWayDopplerObservationSettings->uplinkOneWayDopplerSettings_;
                jsonObject[ K::downlinkOneWayDopplerSettings ] =
                        twoWayDopplerObservationSettings->downlinkOneWayDopplerSettings_;
            }
            return;
        }
        default: {
            return;
        }
    }
}

void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< ObservationSettings >& observationSettings )
{
    using namespace json_interface;
    using K = Keys::Observation;

    const ObservableType observableType =
            getValue< ObservableType >( jsonObject, K::observableType );

    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList =
            getValue< std::vector< std::shared_ptr< LightTimeCorrectionSettings > > >(
                    jsonObject,
                    K::lightTimeCorrectionSettingsList,
                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) );
    std::shared_ptr< ObservationBiasSettings > biasSettings =
            getValue< std::shared_ptr< ObservationBiasSettings > >(
                    jsonObject, K::biasSettings, nullptr );

    switch( observableType )
    {
        case one_way_doppler: {
            std::shared_ptr< DopplerProperTimeRateSettings > transmitterProperTimeRateSettings =
                    getValue< std::shared_ptr< DopplerProperTimeRateSettings > >(
                            jsonObject, K::transmitterProperTimeRateSettings, nullptr );
            std::shared_ptr< DopplerProperTimeRateSettings > receiverProperTimeRateSettings =
                    getValue< std::shared_ptr< DopplerProperTimeRateSettings > >(
                            jsonObject, K::receiverProperTimeRateSettings, nullptr );

            observationSettings = std::make_shared< OneWayDopplerObservationSettings >(
                    lightTimeCorrectionsList,
                    transmitterProperTimeRateSettings,
                    receiverProperTimeRateSettings,
                    biasSettings );
            return;
        }
        case one_way_differenced_range: {
            double constantIntegrationTime =
                    getValue< double >( jsonObject, K::constantIntegrationTime );

            observationSettings = std::make_shared< OneWayDifferencedRangeRateObservationSettings >(
                    [ = ]( const double ) { return constantIntegrationTime; },
                    lightTimeCorrectionsList,
                    biasSettings );
            return;
        }
        case n_way_range: {
            std::vector< std::shared_ptr< ObservationSettings > > oneWayRangeObservationSettings =
                    getValue< std::vector< std::shared_ptr< ObservationSettings > > >(
                            jsonObject,
                            K::oneWayRangeObservationSettings,
                            std::vector< std::shared_ptr< ObservationSettings > >( ) );
            if( oneWayRangeObservationSettings.size( ) == 0 )
            {
                observationSettings = std::make_shared< ObservationSettings >(
                        observableType, lightTimeCorrectionsList, biasSettings );
            }
            else
            {
                std::vector< double > retransmissionTimes = getValue< std::vector< double > >(
                        jsonObject, K::retransmissionTimes, std::vector< double >( ) );
                observationSettings = std::make_shared< NWayRangeObservationSettings >(
                        oneWayRangeObservationSettings,
                        [ = ]( const double ) { return retransmissionTimes; },
                        biasSettings );
            }
            return;
        }
        case two_way_doppler: {
            std::shared_ptr< ObservationSettings > uplinkOneWayDopplerSettings =
                    getValue< std::shared_ptr< ObservationSettings > >(
                            jsonObject, K::uplinkOneWayDopplerSettings, nullptr );
            std::shared_ptr< ObservationSettings > downlinkOneWayDopplerSettings =
                    getValue< std::shared_ptr< ObservationSettings > >(
                            jsonObject, K::downlinkOneWayDopplerSettings, nullptr );
            if( ( uplinkOneWayDopplerSettings == nullptr ) &&
                ( downlinkOneWayDopplerSettings == nullptr ) )
            {
                observationSettings = std::make_shared< ObservationSettings >(
                        observableType, lightTimeCorrectionsList, biasSettings );
            }
            else if( ( uplinkOneWayDopplerSettings != nullptr ) &&
                     ( downlinkOneWayDopplerSettings != nullptr ) )
            {
                if( std::dynamic_pointer_cast< OneWayDopplerObservationSettings >(
                            uplinkOneWayDopplerSettings ) == nullptr ||
                    std::dynamic_pointer_cast< OneWayDopplerObservationSettings >(
                            downlinkOneWayDopplerSettings ) == nullptr )
                {
                    throw std::runtime_error(
                            "Error when reading 2-way Doppler dettings from JSON file, link "
                            "doppler settings are inconsistent" );
                }
                observationSettings = std::make_shared< TwoWayDopplerObservationSettings >(
                        std::dynamic_pointer_cast< OneWayDopplerObservationSettings >(
                                uplinkOneWayDopplerSettings ),
                        std::dynamic_pointer_cast< OneWayDopplerObservationSettings >(
                                downlinkOneWayDopplerSettings ),
                        biasSettings );
            }
            else
            {
                throw std::runtime_error(
                        "Error when reading 2-way Doppler dettings from JSON file, settings are "
                        "inconsistent" );
            }
            return;
        }
        default: {
            observationSettings = std::make_shared< ObservationSettings >(
                    observableType, lightTimeCorrectionsList, biasSettings );
            return;
        }
    }
}

//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< ObservationBiasSettings >& biasSettings )
{
    if( !biasSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::ObservationBias;

    const ObservationBiasTypes biasType = biasSettings->observationBiasType_;
    ;
    jsonObject[ K::biasType ] = biasType;

    switch( biasType )
    {
        case multiple_observation_biases: {
            std::shared_ptr< MultipleObservationBiasSettings > multipleBiasSettings =
                    std::dynamic_pointer_cast< MultipleObservationBiasSettings >( biasSettings );
            assertNonnullptrPointer( multipleBiasSettings );
            jsonObject[ K::multipleBiasesList ] = multipleBiasSettings->biasSettingsList_;
            return;
        }
        case constant_absolute_bias: {
            std::shared_ptr< ConstantObservationBiasSettings > constantBiasSettings =
                    std::dynamic_pointer_cast< ConstantObservationBiasSettings >( biasSettings );
            assertNonnullptrPointer( constantBiasSettings );

            jsonObject[ K::constantBias ] = constantBiasSettings->observationBias_;
            return;
        }
        case constant_relative_bias: {
            std::shared_ptr< ConstantObservationBiasSettings > constantBiasSettings =
                    std::dynamic_pointer_cast< ConstantObservationBiasSettings >( biasSettings );
            assertNonnullptrPointer( constantBiasSettings );

            jsonObject[ K::constantBias ] = constantBiasSettings->observationBias_;
            return;
        }
        case arc_wise_constant_absolute_bias: {
            std::shared_ptr< ArcWiseConstantObservationBiasSettings > arcWiseBiasSettings =
                    std::dynamic_pointer_cast< ArcWiseConstantObservationBiasSettings >(
                            biasSettings );
            assertNonnullptrPointer( arcWiseBiasSettings );

            jsonObject[ K::arcWiseBiasList ] = arcWiseBiasSettings->observationBiases_;
            jsonObject[ K::arcStartTimes ] = arcWiseBiasSettings->arcStartTimes_;
            jsonObject[ K::referenceLinkEnd ] = arcWiseBiasSettings->linkEndForTime_;

            return;
        }
        case arc_wise_constant_relative_bias: {
            std::shared_ptr< ArcWiseConstantObservationBiasSettings > arcWiseBiasSettings =
                    std::dynamic_pointer_cast< ArcWiseConstantObservationBiasSettings >(
                            biasSettings );
            assertNonnullptrPointer( arcWiseBiasSettings );

            jsonObject[ K::arcWiseBiasList ] = arcWiseBiasSettings->observationBiases_;
            jsonObject[ K::arcStartTimes ] = arcWiseBiasSettings->arcStartTimes_;
            jsonObject[ K::referenceLinkEnd ] = arcWiseBiasSettings->linkEndForTime_;

            return;
        }
        default:
            return;
    }
}

void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< ObservationBiasSettings >& biasSettings )
{
    using namespace json_interface;
    using K = Keys::ObservationBias;

    const ObservationBiasTypes biasType =
            getValue< ObservationBiasTypes >( jsonObject, K::biasType );

    switch( biasType )
    {
        case constant_absolute_bias: {
            Eigen::VectorXd constantBias =
                    getValue< Eigen::VectorXd >( jsonObject, K::constantBias );
            biasSettings =
                    std::make_shared< ConstantObservationBiasSettings >( constantBias, true );
            return;
        }
        case constant_relative_bias: {
            Eigen::VectorXd constantBias =
                    getValue< Eigen::VectorXd >( jsonObject, K::constantBias );
            biasSettings =
                    std::make_shared< ConstantObservationBiasSettings >( constantBias, false );
            return;
        }
        case arc_wise_constant_absolute_bias: {
            std::vector< double > arcStartTimes =
                    getValue< std::vector< double > >( jsonObject, K::arcStartTimes );
            std::vector< Eigen::VectorXd > observationBiases =
                    getValue< std::vector< Eigen::VectorXd > >( jsonObject, K::arcWiseBiasList );
            LinkEndType linkEndForTime = getValue< LinkEndType >( jsonObject, K::referenceLinkEnd );
            biasSettings = std::make_shared< ArcWiseConstantObservationBiasSettings >(
                    arcStartTimes, observationBiases, linkEndForTime, true );
            return;
        }
        case arc_wise_constant_relative_bias: {
            std::vector< double > arcStartTimes =
                    getValue< std::vector< double > >( jsonObject, K::arcStartTimes );
            std::vector< Eigen::VectorXd > observationBiases =
                    getValue< std::vector< Eigen::VectorXd > >( jsonObject, K::arcWiseBiasList );
            LinkEndType linkEndForTime = getValue< LinkEndType >( jsonObject, K::referenceLinkEnd );
            biasSettings = std::make_shared< ArcWiseConstantObservationBiasSettings >(
                    arcStartTimes, observationBiases, linkEndForTime, false );
            return;
        }
        case multiple_observation_biases: {
            std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList =
                    getValue< std::vector< std::shared_ptr< ObservationBiasSettings > > >(
                            jsonObject, K::multipleBiasesList );
            biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );

            return;
        }
        default:
            return;
    }
}

//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< DopplerProperTimeRateSettings >& properTimeRateSettings )
{
    if( !properTimeRateSettings )
    {
        return;
    }

    using namespace json_interface;
    using K = Keys::Observation;

    const DopplerProperTimeRateType rateType = properTimeRateSettings->dopplerProperTimeRateType_;
    jsonObject[ K::properTimeRateType ] = rateType;
    switch( rateType )
    {
        case direct_first_order_doppler_proper_time_rate: {
            std::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > firstOrderRateType =
                    std::dynamic_pointer_cast< DirectFirstOrderDopplerProperTimeRateSettings >(
                            properTimeRateSettings );
            assertNonnullptrPointer( firstOrderRateType );
            jsonObject[ K::centralBody ] = firstOrderRateType->centralBodyName_;
            return;
        }
        default:
            return;
    }
}

void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< DopplerProperTimeRateSettings >& properTimeRateSettings )
{
    using namespace json_interface;
    using K = Keys::Observation;

    const DopplerProperTimeRateType rateType =
            getValue< DopplerProperTimeRateType >( jsonObject, K::properTimeRateType );

    switch( rateType )
    {
        case direct_first_order_doppler_proper_time_rate: {
            properTimeRateSettings =
                    std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >(
                            getValue< std::string >( jsonObject, K::centralBody ) );
            return;
        }
        default: {
            return;
        }
    }
}

//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< LightTimeCorrectionSettings >& lightTimeCorrectionSettings )
{
    if( !lightTimeCorrectionSettings )
    {
        return;
    }

    using namespace json_interface;
    using K = Keys::Observation;

    const LightTimeCorrectionType correctionType =
            lightTimeCorrectionSettings->getCorrectionType( );
    jsonObject[ K::lightTimeCorrectionType ] = correctionType;

    switch( correctionType )
    {
        case first_order_relativistic: {
            std::shared_ptr< FirstOrderRelativisticLightTimeCorrectionSettings >
                    firstOrderCorrectionSettings = std::dynamic_pointer_cast<
                            FirstOrderRelativisticLightTimeCorrectionSettings >(
                            lightTimeCorrectionSettings );
            assertNonnullptrPointer( firstOrderCorrectionSettings );
            jsonObject[ K::perturbingBodies ] =
                    firstOrderCorrectionSettings->getPerturbingBodies( );

            return;
        }
        default:
            return;
    }
}

void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< LightTimeCorrectionSettings >& lightTimeCorrectionSettings )
{
    using namespace json_interface;
    using K = Keys::Observation;

    try
    {
        const LightTimeCorrectionType correctionType =
                getValue< LightTimeCorrectionType >( jsonObject, K::lightTimeCorrectionType );

        switch( correctionType )
        {
            case first_order_relativistic: {
                lightTimeCorrectionSettings =
                        std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                getValue< std::vector< std::string > >( jsonObject,
                                                                        K::perturbingBodies ) );
                return;
            }
            default: {
                return;
            }
        }
    }
    catch( ... )
    { }
}

//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< ObservationSimulationTimeSettings< double > >&
                      observationSimulationTimeSettings )
{
    if( !observationSimulationTimeSettings )
    {
        return;
    }

    using namespace json_interface;
    using K = Keys::Observation;

    if( std::dynamic_pointer_cast< TabulatedObservationSimulationTimeSettings< double > >(
                observationSimulationTimeSettings ) != nullptr )
    {
        std::vector< double > simulationTimes =
                std::dynamic_pointer_cast< TabulatedObservationSimulationTimeSettings< double > >(
                        observationSimulationTimeSettings )
                        ->simulationTimes_;
        jsonObject[ K::observationSimulationTimesType ] = tabulated_observation_simulation_times;
        jsonObject[ K::observationSimulationTimesList ] = simulationTimes;

        return;
    }
    else
    {
        return;
    }
}

//! Create a `json` object from a shared pointer to a `ObservationViabilitySettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< ObservationViabilitySettings >& observationViabilitySettings )
{
    if( !observationViabilitySettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Observation;

    const ObservationViabilityType viabilityType =
            observationViabilitySettings->observationViabilityType_;
    jsonObject[ K::observableViabilityType ] = viabilityType;
    jsonObject[ K::associatedLinkEnd ] = observationViabilitySettings->getAssociatedLinkEnd( );
    if( observationViabilitySettings->getDoubleParameter( ) ==
        observationViabilitySettings->getDoubleParameter( ) )
    {
        jsonObject[ K::doubleParameter ] = observationViabilitySettings->getDoubleParameter( );
    }
    if( observationViabilitySettings->getStringParameter( ).size( ) > 0 )
    {
        jsonObject[ K::stringParameter ] = observationViabilitySettings->getStringParameter( );
    }
}

}  // namespace observation_models

}  // namespace tudat

namespace boost
{

template<>
tudat::observation_models::ObservableType lexical_cast( const std::string& s )
{
    return tudat::observation_models::observationTypesInverse.at( s );
}

template<>
std::string lexical_cast( const tudat::observation_models::ObservableType& component )
{
    return tudat::observation_models::observationTypes.at( component );
}

template<>
tudat::observation_models::LinkEndType lexical_cast( const std::string& s )
{
    return tudat::observation_models::linkEndTypesInverse.at( s );
}

template<>
std::string lexical_cast( const tudat::observation_models::LinkEndType& component )
{
    return tudat::observation_models::linkEndTypes.at( component );
}

template<>
tudat::observation_models::LinkEnds lexical_cast( const std::string& s )
{
    tudat::observation_models::LinkEnds linkEnds;
    std::vector< std::string > individualLinkEnds = tudat::json_interface::split( s, ';' );
    for( unsigned int i = 0; i < individualLinkEnds.size( ); i++ )
    {
        std::vector< std::string > typeAndIdentifier =
                tudat::json_interface::split( individualLinkEnds.at( i ), ':' );

        std::vector< std::string > splitLinkEnd = tudat::json_interface::split(
                typeAndIdentifier.at( 1 ).substr( 1, typeAndIdentifier.at( 1 ).size( ) - 2 ),
                ',',
                false );
        if( splitLinkEnd.size( ) == 1 )
        {
            splitLinkEnd.push_back( "" );
        }

        linkEnds[ tudat::observation_models::linkEndTypesInverse.at( typeAndIdentifier.at( 0 ) ) ] =
                std::make_pair( splitLinkEnd.at( 0 ), splitLinkEnd.at( 1 ) );
    }

    return linkEnds;
}

template<>
std::string lexical_cast( const tudat::observation_models::LinkEnds& component )
{
    std::string linkEndsString;

    for( auto it = component.begin( ); it != component.end( ); it++ )
    {
        if( it != component.begin( ) )
        {
            linkEndsString += "; ";
        }

        linkEndsString += boost::lexical_cast< std::string >( it->first ) + ":(" +
                it->second.first + "," + it->second.second + ")";
    }
    return linkEndsString;
}

}  // namespace boost
