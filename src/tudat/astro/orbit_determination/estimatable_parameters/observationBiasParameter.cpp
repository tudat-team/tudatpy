/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/estimatable_parameters/observationBiasParameter.h"

#include <stdexcept>

namespace tudat
{

namespace estimatable_parameters
{

std::string createObsBiasSecondaryIdentifier( const observation_models::ObservableType observableType,
                                              const observation_models::LinkEnds& linkEnds )
{
    std::string transmitterStr = "None";
    std::string receiverStr = "None";

    auto transmitterIt = linkEnds.find( observation_models::transmitter );
    if( transmitterIt != linkEnds.end( ) )
    {
        transmitterStr = transmitterIt->second.getReferencePointName( );
    }

    auto receiverIt = linkEnds.find( observation_models::receiver );
    if( receiverIt != linkEnds.end( ) )
    {
        receiverStr = receiverIt->second.getReferencePointName( );
    }

    return transmitterStr + " --> " + receiverStr + " , " + observation_models::getObservableName( observableType );
}

SingleArcObservationBiasParameter::SingleArcObservationBiasParameter(
        const EstimatebleParametersEnum parameterName,
        const std::function< Eigen::VectorXd( ) > getCurrentBias,
        const std::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::string& pointOnBodyId ):
    ObservationBiasFunctionWrapper< Eigen::VectorXd >(
            parameterName,
            linkEnds.begin( )->second.bodyName_,
            pointOnBodyId.empty( ) ? createObsBiasSecondaryIdentifier( observableType, linkEnds ) : pointOnBodyId,
            getCurrentBias,
            resetCurrentBias ),
    linkEnds_( linkEnds ),
    observableType_( observableType )
{ }

Eigen::VectorXd SingleArcObservationBiasParameter::getParameterValue( )
{
    if( biasFunctionsAreDefined( ) )
    {
        return getBiasFunction( )( );
    }
    else if( hasDeferredBiasValue( ) )
    {
        return getDeferredBiasValue( );
    }
    else
    {
        return Eigen::VectorXd::Constant( getParameterSize( ), TUDAT_NAN );
    }
}

void SingleArcObservationBiasParameter::setParameterValue( Eigen::VectorXd parameterValue )
{
    if( getParameterSize( ) != parameterValue.rows( ) )
    {
        throw std::runtime_error( "Error, size of parameter incompatible with expected size when resetting value." );
    }

    resetOrDeferBiasValue( parameterValue );
}

int SingleArcObservationBiasParameter::getParameterSize( )
{
    return observation_models::getObservableSize( observableType_ );
}

void SingleArcObservationBiasParameter::throwExceptionIfNotFullyDefined( )
{
    if( !biasFunctionsAreDefined( ) )
    {
        throw std::runtime_error(
                "Error in " + getParameterTypeString( parameterName_.first ) + " of observable type " +
                observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + " with link ends: " +
                observation_models::getLinkEndsString( linkEnds_ ) +
                " parameter not linked to bias object. Associated bias model been implemented in observation model. "
                "This may be because you are resetting the parameter value before creating observation models, "
                "or because you have not defined the required bias model." );
    }
}

observation_models::LinkEnds SingleArcObservationBiasParameter::getLinkEnds( )
{
    return linkEnds_;
}

observation_models::ObservableType SingleArcObservationBiasParameter::getObservableType( )
{
    return observableType_;
}

std::string SingleArcObservationBiasParameter::getParameterDescription( )
{
    std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
            observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + ") and link ends: (" +
            observation_models::getLinkEndsString( linkEnds_ ) + ")";
    return parameterDescription;
}

MultiArcObservationBiasParameter::MultiArcObservationBiasParameter(
        const EstimatebleParametersEnum parameterName,
        const std::vector< double > arcStartTimes,
        const std::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
        const std::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
        const int linkEndIndex,
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::string& pointOnBodyId ):
    ObservationBiasFunctionWrapper< std::vector< Eigen::VectorXd > >(
            parameterName,
            linkEnds.begin( )->second.bodyName_,
            pointOnBodyId.empty( ) ? createObsBiasSecondaryIdentifier( observableType, linkEnds ) : pointOnBodyId,
            getBiasList,
            resetBiasList ),
    arcStartTimes_( arcStartTimes ),
    linkEndIndex_( linkEndIndex ),
    linkEnds_( linkEnds ),
    observableType_( observableType ),
    observableSize_( observation_models::getObservableSize( observableType ) ),
    numberOfArcs_( static_cast< int >( arcStartTimes.size( ) ) )
{ }

Eigen::VectorXd MultiArcObservationBiasParameter::getParameterValue( )
{
    if( biasFunctionsAreDefined( ) )
    {
        std::vector< Eigen::VectorXd > observationBiases = getBiasFunction( )( );
        Eigen::VectorXd currentParameterSet = Eigen::VectorXd::Zero( observableSize_ * observationBiases.size( ) );
        for( unsigned int i = 0; i < observationBiases.size( ); i++ )
        {
            currentParameterSet.segment( i * observableSize_, observableSize_ ) = observationBiases.at( i );
        }
        return currentParameterSet;
    }
    else if( hasDeferredBiasValue( ) )
    {
        const std::vector< Eigen::VectorXd >& observationBiases = getDeferredBiasValue( );
        Eigen::VectorXd currentParameterSet = Eigen::VectorXd::Zero( observableSize_ * observationBiases.size( ) );
        for( unsigned int i = 0; i < observationBiases.size( ); i++ )
        {
            currentParameterSet.segment( i * observableSize_, observableSize_ ) = observationBiases.at( i );
        }
        return currentParameterSet;
    }
    else
    {
        return Eigen::VectorXd::Constant( getParameterSize( ), TUDAT_NAN );
    }
}

void MultiArcObservationBiasParameter::setParameterValue( Eigen::VectorXd parameterValue )
{
    if( getParameterSize( ) != parameterValue.rows( ) )
    {
        throw std::runtime_error( "Error, size of parameter incompatible with expected size when resetting value." );
    }

    std::vector< Eigen::VectorXd > observationBiases;
    for( int i = 0; i < numberOfArcs_; i++ )
    {
        observationBiases.push_back( parameterValue.segment( i * observableSize_, observableSize_ ) );
    }

    resetOrDeferBiasValue( observationBiases );
}

int MultiArcObservationBiasParameter::getParameterSize( )
{
    return observableSize_ * numberOfArcs_;
}

void MultiArcObservationBiasParameter::throwExceptionIfNotFullyDefined( )
{
    if( !biasFunctionsAreDefined( ) )
    {
        throw std::runtime_error(
                "Error in " + getParameterTypeString( parameterName_.first ) + " of observable type " +
                observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + " with link ends: " +
                observation_models::getLinkEndsString( linkEnds_ ) +
                " parameter not linked to bias object. Has associated bias model been implemented in observation model?" );
    }
}

observation_models::LinkEnds MultiArcObservationBiasParameter::getLinkEnds( )
{
    return linkEnds_;
}

observation_models::ObservableType MultiArcObservationBiasParameter::getObservableType( )
{
    return observableType_;
}

std::string MultiArcObservationBiasParameter::getParameterDescription( )
{
    std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
            observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + ") and link ends: (" +
            observation_models::getLinkEndsString( linkEnds_ ) + ")";
    return parameterDescription;
}

std::vector< double > MultiArcObservationBiasParameter::getArcStartTimes( )
{
    return arcStartTimes_;
}

int MultiArcObservationBiasParameter::getLinkEndIndex( )
{
    return linkEndIndex_;
}

std::shared_ptr< interpolators::LookUpScheme< double > > MultiArcObservationBiasParameter::getLookupScheme( )
{
    return lookupScheme_;
}

void MultiArcObservationBiasParameter::setLookupScheme(
        const std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme )
{
    lookupScheme_ = lookupScheme;
}

void TimeBiasParameterBase::setBodyAccelerationFunction(
        const std::function< Eigen::VectorXd( const double ) > bodyAccelerationFunction )
{
    bodyAccelerationFunction_ = bodyAccelerationFunction;
}

std::function< Eigen::VectorXd( const double ) > TimeBiasParameterBase::getBodyAccelerationFunction( )
{
    return bodyAccelerationFunction_;
}

SingleArcTimeBiasParameter::SingleArcTimeBiasParameter(
        const EstimatebleParametersEnum parameterName,
        const std::function< Eigen::VectorXd( ) > getCurrentBias,
        const std::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
        const observation_models::LinkEndType linkEndForTime,
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::string& pointOnBodyId ):
    ObservationBiasFunctionWrapper< Eigen::VectorXd >(
            parameterName,
            linkEnds.begin( )->second.bodyName_,
            pointOnBodyId.empty( ) ? createObsBiasSecondaryIdentifier( observableType, linkEnds ) : pointOnBodyId,
            getCurrentBias,
            resetCurrentBias ),
    linkEndForTime_( linkEndForTime ),
    linkEnds_( linkEnds ),
    observableType_( observableType )
{
    linkEndIndex_ = observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                            observableType_, linkEndForTime_, linkEnds_.size( ) )
                            .at( 0 );
}

Eigen::VectorXd SingleArcTimeBiasParameter::getParameterValue( )
{
    if( biasFunctionsAreDefined( ) )
    {
        return getBiasFunction( )( );
    }
    else if( hasDeferredBiasValue( ) )
    {
        return getDeferredBiasValue( );
    }
    else
    {
        return Eigen::VectorXd::Constant( 1, TUDAT_NAN );
    }
}

void SingleArcTimeBiasParameter::setParameterValue( Eigen::VectorXd parameterValue )
{
    if( getParameterSize( ) != parameterValue.rows( ) )
    {
        throw std::runtime_error(
                "Error, size of parameter (type:constant_time_observation_bias) incompatible with expected size when resetting "
                "value." );
    }

    resetOrDeferBiasValue( parameterValue );
}

int SingleArcTimeBiasParameter::getParameterSize( )
{
    return 1;
}

void SingleArcTimeBiasParameter::throwExceptionIfNotFullyDefined( )
{
    if( !biasFunctionsAreDefined( ) )
    {
        throw std::runtime_error(
                "Error in " + getParameterTypeString( parameterName_.first ) + " of observable type " +
                observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + " with link ends: " +
                observation_models::getLinkEndsString( linkEnds_ ) +
                " parameter not linked to bias object. Associated bias model been implemented in observation model. "
                "This may be because you are resetting the parameter value before creating observation models, "
                "or because you have not defined the required bias model." );
    }
}

observation_models::LinkEnds SingleArcTimeBiasParameter::getLinkEnds( )
{
    return linkEnds_;
}

observation_models::LinkEndId SingleArcTimeBiasParameter::getLinkEndId( )
{
    return linkEnds_.at( linkEndForTime_ );
}

observation_models::LinkEndType SingleArcTimeBiasParameter::getReferenceLinkEnd( )
{
    return linkEndForTime_;
}

observation_models::ObservableType SingleArcTimeBiasParameter::getObservableType( )
{
    return observableType_;
}

std::string SingleArcTimeBiasParameter::getParameterDescription( )
{
    std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
            observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + ") and link ends: (" +
            observation_models::getLinkEndsString( linkEnds_ ) + ")";
    return parameterDescription;
}

int SingleArcTimeBiasParameter::getLinkEndIndex( )
{
    return linkEndIndex_;
}

MultiArcTimeBiasParameter::MultiArcTimeBiasParameter(
        const EstimatebleParametersEnum parameterName,
        const std::vector< double > arcStartTimes,
        const std::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
        const std::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
        const observation_models::LinkEndType linkEndForTime,
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::string& pointOnBodyId ):
    ObservationBiasFunctionWrapper< std::vector< Eigen::VectorXd > >(
            parameterName,
            linkEnds.begin( )->second.bodyName_,
            pointOnBodyId.empty( ) ? createObsBiasSecondaryIdentifier( observableType, linkEnds ) : pointOnBodyId,
            getBiasList,
            resetBiasList ),
    arcStartTimes_( arcStartTimes ),
    linkEndForTime_( linkEndForTime ),
    linkEnds_( linkEnds ),
    observableType_( observableType ),
    numberOfArcs_( static_cast< int >( arcStartTimes.size( ) ) )
{
    linkEndIndex_ = observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                            observableType_, linkEndForTime_, linkEnds_.size( ) )
                            .at( 0 );
}

Eigen::VectorXd MultiArcTimeBiasParameter::getParameterValue( )
{
    if( biasFunctionsAreDefined( ) )
    {
        std::vector< Eigen::VectorXd > observationBiases = getBiasFunction( )( );
        Eigen::VectorXd currentParameterSet = Eigen::VectorXd::Zero( observationBiases.size( ) );
        for( unsigned int i = 0; i < observationBiases.size( ); i++ )
        {
            currentParameterSet.segment( i, 1 ) = observationBiases.at( i );
        }
        return currentParameterSet;
    }
    else if( hasDeferredBiasValue( ) )
    {
        const std::vector< Eigen::VectorXd >& observationBiases = getDeferredBiasValue( );
        Eigen::VectorXd currentParameterSet = Eigen::VectorXd::Zero( observationBiases.size( ) );
        for( unsigned int i = 0; i < observationBiases.size( ); i++ )
        {
            currentParameterSet.segment( i, 1 ) = observationBiases.at( i );
        }
        return currentParameterSet;
    }
    else
    {
        return Eigen::VectorXd::Constant( getParameterSize( ), TUDAT_NAN );
    }
}

void MultiArcTimeBiasParameter::setParameterValue( Eigen::VectorXd parameterValue )
{
    if( getParameterSize( ) != parameterValue.rows( ) )
    {
        throw std::runtime_error(
                "Error, size of parameter (type:arc_wise_time_observation_bias) incompatible with expected size when resetting "
                "value." );
    }

    std::vector< Eigen::VectorXd > observationBiases;
    for( int i = 0; i < numberOfArcs_; i++ )
    {
        observationBiases.push_back( parameterValue.segment( i, 1 ) );
    }

    resetOrDeferBiasValue( observationBiases );
}

int MultiArcTimeBiasParameter::getParameterSize( )
{
    return numberOfArcs_;
}

void MultiArcTimeBiasParameter::throwExceptionIfNotFullyDefined( )
{
    if( !biasFunctionsAreDefined( ) )
    {
        throw std::runtime_error(
                "Error in " + getParameterTypeString( parameterName_.first ) + " of observable type " +
                observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + " with link ends: " +
                observation_models::getLinkEndsString( linkEnds_ ) +
                " parameter not linked to bias object. Associated bias model been implemented in observation model. "
                "This may be because you are resetting the parameter value before creating observation models, "
                "or because you have not defined the required bias model." );
    }
}

observation_models::LinkEnds MultiArcTimeBiasParameter::getLinkEnds( )
{
    return linkEnds_;
}

observation_models::ObservableType MultiArcTimeBiasParameter::getObservableType( )
{
    return observableType_;
}

std::string MultiArcTimeBiasParameter::getParameterDescription( )
{
    std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
            observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + ") and link ends: (" +
            observation_models::getLinkEndsString( linkEnds_ ) + ")";
    return parameterDescription;
}

std::vector< double > MultiArcTimeBiasParameter::getArcStartTimes( )
{
    return arcStartTimes_;
}

int MultiArcTimeBiasParameter::getLinkEndIndex( )
{
    return linkEndIndex_;
}

std::shared_ptr< interpolators::LookUpScheme< double > > MultiArcTimeBiasParameter::getLookupScheme( )
{
    return lookupScheme_;
}

void MultiArcTimeBiasParameter::setLookupScheme( const std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme )
{
    lookupScheme_ = lookupScheme;
}

observation_models::LinkEndId MultiArcTimeBiasParameter::getLinkEndId( )
{
    return linkEnds_.at( linkEndForTime_ );
}

observation_models::LinkEndType MultiArcTimeBiasParameter::getReferenceLinkEnd( )
{
    return linkEndForTime_;
}

ConstantTimeDriftBiasParameter::ConstantTimeDriftBiasParameter(
        const EstimatebleParametersEnum parameterName,
        const std::function< Eigen::VectorXd( ) > getCurrentBias,
        const std::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
        const int linkEndIndex,
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const double referenceEpoch,
        const std::string& pointOnBodyId ):
    SingleArcObservationBiasParameter(
            parameterName, getCurrentBias, resetCurrentBias, linkEnds, observableType, pointOnBodyId ),
    linkEndIndex_( linkEndIndex ),
    referenceEpoch_( referenceEpoch )
{ }

int ConstantTimeDriftBiasParameter::getLinkEndIndex( )
{
    return linkEndIndex_;
}

double ConstantTimeDriftBiasParameter::getReferenceEpoch( )
{
    return referenceEpoch_;
}

ArcWiseTimeDriftBiasParameter::ArcWiseTimeDriftBiasParameter(
        const EstimatebleParametersEnum parameterName,
        const std::vector< double > arcStartTimes,
        const std::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
        const std::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
        const int linkEndIndex,
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const std::vector< double > referenceEpochs,
        const std::string& pointOnBodyId ):
    MultiArcObservationBiasParameter( parameterName,
                                      arcStartTimes,
                                      getBiasList,
                                      resetBiasList,
                                      linkEndIndex,
                                      linkEnds,
                                      observableType,
                                      pointOnBodyId ),
    referenceEpochs_( referenceEpochs )
{ }

std::vector< double > ArcWiseTimeDriftBiasParameter::getReferenceEpochs( )
{
    return referenceEpochs_;
}

}  // namespace estimatable_parameters

}  // namespace tudat
