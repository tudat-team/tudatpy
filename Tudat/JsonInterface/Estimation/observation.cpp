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
#include "Tudat/JsonInterface/Support/utilities.h"

namespace tudat
{

namespace observation_models
{
// SingleDependentVariableSaveSettings

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( nlohmann::json& jsonObject,
              const boost::shared_ptr< ObservationSettings >& observationSettings )
{
    if ( !observationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Observation;

    const ObservableType observableType = observationSettings->observableType_;
    jsonObject[ K::observableType ] = observableType;
    assignIfNotEmpty( jsonObject, K::lightTimeCorrectionSettingsList, observationSettings->lightTimeCorrectionsList_ );
    if( observationSettings->biasSettings_ != NULL )
    {
        jsonObject[ K::biasSettings ] = observationSettings->biasSettings_;
    }

    switch ( observableType )
    {
    case one_way_range:
    case angular_position:
    case position_observable:
    case one_way_doppler:
    {
        boost::shared_ptr< OneWayDopplerObservationSettings > oneWayDopplerObservationSettings =
                boost::dynamic_pointer_cast< OneWayDopplerObservationSettings >( observationSettings );
        if( !( oneWayDopplerObservationSettings == NULL ) )
        {
            if( oneWayDopplerObservationSettings->transmitterProperTimeRateSettings_ != NULL )
            {
                jsonObject[ K::transmitterProperTimeRateSettings ] = oneWayDopplerObservationSettings->transmitterProperTimeRateSettings_;
            }

            if( oneWayDopplerObservationSettings->receiverProperTimeRateSettings_ != NULL )
            {
                jsonObject[ K::receiverProperTimeRateSettings ] = oneWayDopplerObservationSettings->receiverProperTimeRateSettings_ ;
            }
            return;
        }
        return;
    }
    case one_way_differenced_range:
    {
        boost::shared_ptr< OneWayDifferencedRangeRateObservationSettings > oneWayDifferencedRangeObservationSettings =
                boost::dynamic_pointer_cast< OneWayDifferencedRangeRateObservationSettings >( observationSettings );
        assertNonNullPointer( oneWayDifferencedRangeObservationSettings );

        std::cerr<<"Warning, assuming constant integration time of 1-way differenced range when parsing to JSON, correctness not ensured"<<std::endl;
        jsonObject[ K::constantIntegrationTime ] = oneWayDifferencedRangeObservationSettings->integrationTimeFunction_( 0.0 );
        return;
    }
    case n_way_range:
    {
        boost::shared_ptr< NWayRangeObservationSettings > nWayRangeObservationSettings =
                boost::dynamic_pointer_cast< NWayRangeObservationSettings >( observationSettings );
        if( !( nWayRangeObservationSettings == NULL ) )
        {
            jsonObject[ K::oneWayRangeObsevationSettings ] = nWayRangeObservationSettings->oneWayRangeObsevationSettings_;

            if( !nWayRangeObservationSettings->retransmissionTimesFunction_.empty( ) )
            {
                std::vector< double > retransmissionTimes = nWayRangeObservationSettings->retransmissionTimesFunction_( 0.0 );
                std::cerr<<"Warning, assuming constant delay times of n-way range when parsing to JSON, correctness not ensured"<<std::endl;
                jsonObject[ K::retransmissionTimes ] = retransmissionTimes;
            }
            return;
        }
    }
    case two_way_doppler:
    {
        boost::shared_ptr< TwoWayDopplerObservationSettings > twoWayDopplerObservationSettings =
                boost::dynamic_pointer_cast< TwoWayDopplerObservationSettings >( observationSettings );
        if( !( twoWayDopplerObservationSettings == NULL ) )
        {
            jsonObject[ K::uplinkOneWayDopplerSettings ] = twoWayDopplerObservationSettings->uplinkOneWayDopplerSettings_;
            jsonObject[ K::downlinkOneWayDopplerSettings ] = twoWayDopplerObservationSettings->downlinkOneWayDopplerSettings_;
            return;
        }
    }
    default:
    {
        return;
    }
    }
}

//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ObservationBiasSettings >& biasSettings )
{
    if ( !biasSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::ObservationBias;

    const ObservationBiasTypes biasType = biasSettings->observationBiasType_;;
    jsonObject[ K::biasType ] = biasType;

    switch( biasType )
    {
    case multiple_observation_biases:
    {
        boost::shared_ptr< MultipleObservationBiasSettings > multipleBiasSettings =
                boost::dynamic_pointer_cast< MultipleObservationBiasSettings >( biasSettings );
        assertNonNullPointer( multipleBiasSettings );
        jsonObject[ K::multipleBiasesList ] = multipleBiasSettings->biasSettingsList_;
        return;
    }
    case constant_absolute_bias:
    {
        boost::shared_ptr< ConstantObservationBiasSettings > constantBiasSettings =
                boost::dynamic_pointer_cast< ConstantObservationBiasSettings >( biasSettings );
        assertNonNullPointer( constantBiasSettings );

        jsonObject[ K::constantBias ] = constantBiasSettings->observationBias_;
        return;
    }
    case constant_relative_bias:
    {
        boost::shared_ptr< ConstantObservationBiasSettings > constantBiasSettings =
                boost::dynamic_pointer_cast< ConstantObservationBiasSettings >( biasSettings );
        assertNonNullPointer( constantBiasSettings );

        jsonObject[ K::constantBias ] = constantBiasSettings->observationBias_;
        return;
    }
    case arc_wise_constant_absolute_bias:
    {
        boost::shared_ptr< ArcWiseConstantObservationBiasSettings > arcWiseBiasSettings =
                boost::dynamic_pointer_cast< ArcWiseConstantObservationBiasSettings >( biasSettings );
        assertNonNullPointer( arcWiseBiasSettings );

        jsonObject[ K::arcWiseBiasList ] = arcWiseBiasSettings->observationBiases_;
        jsonObject[ K::arcStartTimes ] = arcWiseBiasSettings->arcStartTimes_;
        jsonObject[ K::referenceLinkEnd ] = arcWiseBiasSettings->linkEndForTime_;

        return;
    }
    case arc_wise_constant_relative_bias:
    {
        boost::shared_ptr< ArcWiseConstantObservationBiasSettings > arcWiseBiasSettings =
                boost::dynamic_pointer_cast< ArcWiseConstantObservationBiasSettings >( biasSettings );
        assertNonNullPointer( arcWiseBiasSettings );

        jsonObject[ K::arcWiseBiasList ] = arcWiseBiasSettings->observationBiases_;
        jsonObject[ K::arcStartTimes ] = arcWiseBiasSettings->arcStartTimes_;
        jsonObject[ K::referenceLinkEnd ] = arcWiseBiasSettings->linkEndForTime_;

        return;
    }
    default:
        return;
    }
}

//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< DopplerProperTimeRateSettings >& properTimeRateSettings )
{
    if ( !properTimeRateSettings )
    {
        return;
    }

    using namespace json_interface;
    using K = Keys::Observation;

    const DopplerProperTimeRateType rateType = properTimeRateSettings->dopplerProperTimeRateType_;
    jsonObject[ K::properTimeRateType ] = rateType;
    switch( rateType )
    {
    case direct_first_order_doppler_proper_time_rate:
    {
        boost::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > firstOrderRateType =
                boost::dynamic_pointer_cast< DirectFirstOrderDopplerProperTimeRateSettings >( properTimeRateSettings );
        assertNonNullPointer( firstOrderRateType );
        jsonObject[ K::centralBody ] = firstOrderRateType->centralBodyName_;
        return;
    }
    default:
        return;
    }
}


//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< LightTimeCorrectionSettings >& lightTimeCorrectionSettings )
{
    if ( !lightTimeCorrectionSettings )
    {
        return;
    }

    using namespace json_interface;
    using K = Keys::Observation;

    const LightTimeCorrectionType correctionType = lightTimeCorrectionSettings->getCorrectionType( );
    jsonObject[ K::lightTimeCorrectionType ] = correctionType;

    switch( correctionType )
    {
    case first_order_relativistic:
    {
        boost::shared_ptr< FirstOrderRelativisticLightTimeCorrectionSettings > firstOrderCorrectionSettings =
                boost::dynamic_pointer_cast< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimeCorrectionSettings );
        assertNonNullPointer( firstOrderCorrectionSettings );
        jsonObject[ K::perturbingBodies ] = firstOrderCorrectionSettings->getPerturbingBodies( );

        return;
    }
    default:
        return;
    }
}

//! Create a `json` object from a shared pointer to a `ObservationSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ObservationSimulationTimeSettings< double > >& observationSimulationTimeSettings )
{
    if ( !observationSimulationTimeSettings )
    {
        return;
    }

    using namespace json_interface;
    using K = Keys::Observation;

    if( boost::dynamic_pointer_cast< TabulatedObservationSimulationTimeSettings< double > >(
                observationSimulationTimeSettings ) != NULL )
    {
        std::vector< double > simulationTimes = boost::dynamic_pointer_cast<
                TabulatedObservationSimulationTimeSettings< double > >( observationSimulationTimeSettings )->simulationTimes_;
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
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ObservationViabilitySettings >& observationViabilitySettings )
{
    if ( !observationViabilitySettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Observation;

    const ObservationViabilityType viabilityType = observationViabilitySettings->observationViabilityType_;
    jsonObject[ K::observableViabilityType ] = viabilityType;
    jsonObject[ K::associatedLinkEnd ] = observationViabilitySettings->getAssociatedLinkEnd( );
    if( observationViabilitySettings->getDoubleParameter( ) == observationViabilitySettings->getDoubleParameter( ) )
    {
        jsonObject[ K::doubleParameter ] = observationViabilitySettings->getDoubleParameter( );
    }
    if( observationViabilitySettings->getStringParameter( ).size( ) > 0 )
    {
        jsonObject[ K::stringParameter ] = observationViabilitySettings->getStringParameter( );
    }
}



}

} // namespace tudat

namespace boost
{

template<>
tudat::observation_models::ObservableType lexical_cast( const std::string & s )
{
    return tudat::observation_models::observationTypesInverse.at( s );
}

template<>
std::string lexical_cast(const tudat::observation_models::ObservableType & component )
{
    return tudat::observation_models::observationTypes.at( component );
}

template<>
tudat::observation_models::LinkEndType lexical_cast( const std::string & s )
{
    return tudat::observation_models::linkEndTypesInverse.at( s );
}


template<>
std::string lexical_cast(const tudat::observation_models::LinkEndType & component )
{
    return tudat::observation_models::linkEndTypes.at( component );
}


template<>
tudat::observation_models::LinkEnds lexical_cast( const std::string & s )
{
    tudat::observation_models::LinkEnds linkEnds;
    std::vector< std::string > individualLinkEnds = tudat::json_interface::split( s, ';' );
    for( unsigned int i = 0; i < individualLinkEnds.size( ); i++ )
    {
        std::vector< std::string > typeAndIdentifier =
                tudat::json_interface::split( individualLinkEnds.at( i ), ':' );

        std::vector< std::string > splitLinkEnd =
                tudat::json_interface::split(
                    typeAndIdentifier.at( 1 ).substr(1, typeAndIdentifier.at( 1 ).size( ) - 2 ), ',', false );

        linkEnds[ tudat::observation_models::linkEndTypesInverse.at( typeAndIdentifier.at( 0 ) ) ] =
                std::make_pair( splitLinkEnd.at( 0 ), splitLinkEnd.at( 1 ) );
    }
    return linkEnds;
}

template<>
std::string lexical_cast(const tudat::observation_models::LinkEnds & component )
{
    std::string linkEndsString;
    for( auto it = component.begin( ); it != component.end( ); it++ )
    {
        if( it != component.begin( ) )
        {
            linkEndsString += "; ";
        }

        linkEndsString += boost::lexical_cast< std::string >( it->first ) + ":(" + it->second.first + ", " + it->second.second + ")";
    }
    return linkEndsString;
}

}
