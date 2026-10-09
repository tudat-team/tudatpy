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
    linkEnds_( linkEnds ), observableType_( observableType )
{}

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
        throw std::runtime_error( "Error in " + getParameterTypeString( parameterName_.first ) + " of observable type " +
                                  observation_models::getObservableName( observableType_, linkEnds_.size( ) ) +
                                  " with link ends: " + observation_models::getLinkEndsString( linkEnds_ ) +
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
    arcStartTimes_( arcStartTimes ), linkEndIndex_( linkEndIndex ), linkEnds_( linkEnds ), observableType_( observableType ),
    observableSize_( observation_models::getObservableSize( observableType ) ), numberOfArcs_( static_cast< int >( arcStartTimes.size( ) ) )
{}

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
                observation_models::getObservableName( observableType_, linkEnds_.size( ) ) +
                " with link ends: " + observation_models::getLinkEndsString( linkEnds_ ) +
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

void MultiArcObservationBiasParameter::setLookupScheme( const std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme )
{
    lookupScheme_ = lookupScheme;
}

//! Define the shared selector and validate the supported bias type and common arc boundaries.
SharedObservationBiasParameter::SharedObservationBiasParameter( const observation_models::ObservationBiasTypes biasType,
                                                                const observation_models::ObservableType observableType,
                                                                const observation_models::LinkEndType linkEndType,
                                                                const observation_models::LinkEndId& linkEndId,
                                                                const std::vector< double >& arcStartTimes,
                                                                const observation_models::LinkEndType timeLinkEnd ):
    EstimatableParameter< Eigen::VectorXd >( shared_observation_bias,
                                             linkEndId.bodyName_,
                                             linkEndId.getReferencePointName( ) + ":" +
                                                     observation_models::getLinkEndTypeString( linkEndType ) + ":" +
                                                     observation_models::getObservableName( observableType ) + ":" +
                                                     observation_models::getObservationBiasTypeString( biasType ) ),
    biasType_( biasType ), observableType_( observableType ), linkEndType_( linkEndType ), linkEndId_( linkEndId ),
    arcStartTimes_( arcStartTimes ), timeLinkEnd_( timeLinkEnd == observation_models::unidentified_link_end ? linkEndType : timeLinkEnd )
{
    observation_models::validateSharedObservationBiasSettings( biasType, linkEndType, linkEndId, arcStartTimes, timeLinkEnd );
}

//! Return one observable-sized bias vector per arc, independent of the number of linked models.
int SharedObservationBiasParameter::getParameterSize( )
{
    return observation_models::getObservableSize( observableType_ ) *
            ( biasType_ == observation_models::arc_wise_constant_absolute_bias ? arcStartTimes_.size( ) : 1 );
}

//! Describe the shared selector and bias type for parameter listings and error messages.
std::string SharedObservationBiasParameter::getParameterDescription( )
{
    return "shared observation bias for " + observation_models::getObservableName( observableType_ ) + ", link-end role " +
            observation_models::getLinkEndTypeString( linkEndType_ ) + ", (" + linkEndId_.bodyName_ + ", " +
            linkEndId_.getReferencePointName( ) + "), bias type " + observation_models::getObservationBiasTypeString( biasType_ );
}

//! Select observations by type and by the complete link-end identifier at the requested role.
bool SharedObservationBiasParameter::doesObservationMatch( const observation_models::LinkEnds& linkEnds,
                                                           const observation_models::ObservableType observableType ) const
{
    const auto it = linkEnds.find( linkEndType_ );
    return observableType == observableType_ && it != linkEnds.end( ) && it->second == linkEndId_;
}

//! Return the shared value while checking that all linked models remain consistent.
Eigen::VectorXd SharedObservationBiasParameter::getParameterValue( )
{
    // Match ordinary bias parameters: values may be assigned before the observation models exist.
    if( members_.empty( ) )
    {
        return hasDeferredValue_ ? deferredValue_ : Eigen::VectorXd::Constant( getParameterSize( ), TUDAT_NAN );
    }
    // A shared parameter represents one value; never silently choose between inconsistent models.
    const Eigen::VectorXd value = members_.begin( )->second->getParameterValue( );
    for( const auto& member : members_ )
    {
        if( !value.isApprox( member.second->getParameterValue( ), 1.0E-14 ) )
        {
            throw std::runtime_error( "Inconsistent values in " + getParameterDescription( ) );
        }
    }
    return value;
}

//! Update bound models, or defer an assignment until the first set of models is linked.
void SharedObservationBiasParameter::setParameterValue( Eigen::VectorXd value )
{
    if( value.size( ) != getParameterSize( ) )
    {
        throw std::runtime_error( "Incorrect shared observation bias parameter size." );
    }
    // Only assignments made before binding initialize new models; estimator updates affect current models only.
    if( members_.empty( ) )
    {
        deferredValue_ = value;
        hasDeferredValue_ = true;
    }
    for( const auto& member : members_ )
    {
        member.second->setParameterValue( value );
    }
}

//! Consume the pre-binding assignment after it has initialized every selected model.
void SharedObservationBiasParameter::completeBinding( )
{
    hasDeferredValue_ = false;
    deferredValue_.resize( 0 );
}

//! Discard old model bindings; subsequent binding reads initial values from the new models.
void SharedObservationBiasParameter::clearMembers( )
{
    members_.clear( );
}

//! Create an unbound ordinary bias parameter for the selected link geometry.
std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > SharedObservationBiasParameter::createMember(
        const observation_models::LinkEnds& linkEnds ) const
{
    using namespace observation_models;
    if( biasType_ == arc_wise_constant_absolute_bias )
    {
        // Convert the common time-link role to the event index for this particular link geometry.
        return std::make_shared< MultiArcObservationBiasParameter >(
                arcwise_constant_additive_observation_bias,
                arcStartTimes_,
                nullptr,
                nullptr,
                getLinkEndIndicesForLinkEndTypeAtObservable( observableType_, timeLinkEnd_, linkEnds.size( ) ).at( 0 ),
                linkEnds,
                observableType_ );
    }
    return std::make_shared< SingleArcObservationBiasParameter >(
            biasType_ == constant_absolute_bias ? constant_additive_observation_bias : constant_relative_observation_bias,
            nullptr,
            nullptr,
            linkEnds,
            observableType_ );
}

//! Register a bound member after checking its size, uniqueness, and initial bias value.
void SharedObservationBiasParameter::addMember( const observation_models::LinkEnds& linkEnds,
                                                const std::shared_ptr< EstimatableParameter< Eigen::VectorXd > >& member )
{
    if( members_.count( linkEnds ) != 0 || member->getParameterSize( ) != getParameterSize( ) )
    {
        throw std::runtime_error( "Duplicate or incompatible member of " + getParameterDescription( ) );
    }
    if( hasDeferredValue_ )
    {
        // An explicit parameter assignment overrides the initial values in the model settings.
        member->setParameterValue( deferredValue_ );
    }
    else if( !members_.empty( ) && !getParameterValue( ).isApprox( member->getParameterValue( ), 1.0E-14 ) )
    {
        throw std::runtime_error( "Shared observation bias models must have equal initial values." );
    }
    members_.emplace( linkEnds, member );
}

//! Look up the ordinary parameter used to create partials for a particular link geometry.
std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > SharedObservationBiasParameter::getMember(
        const observation_models::LinkEnds& linkEnds ) const
{
    const auto it = members_.find( linkEnds );
    return it == members_.end( ) ? nullptr : it->second;
}

void TimeBiasParameterBase::setBodyAccelerationFunction( const std::function< Eigen::VectorXd( const double ) > bodyAccelerationFunction )
{
    bodyAccelerationFunction_ = bodyAccelerationFunction;
}

std::function< Eigen::VectorXd( const double ) > TimeBiasParameterBase::getBodyAccelerationFunction( )
{
    return bodyAccelerationFunction_;
}

SingleArcTimeBiasParameter::SingleArcTimeBiasParameter( const EstimatebleParametersEnum parameterName,
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
    linkEndForTime_( linkEndForTime ), linkEnds_( linkEnds ), observableType_( observableType )
{
    linkEndIndex_ =
            observation_models::getLinkEndIndicesForLinkEndTypeAtObservable( observableType_, linkEndForTime_, linkEnds_.size( ) ).at( 0 );
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
        throw std::runtime_error( "Error in " + getParameterTypeString( parameterName_.first ) + " of observable type " +
                                  observation_models::getObservableName( observableType_, linkEnds_.size( ) ) +
                                  " with link ends: " + observation_models::getLinkEndsString( linkEnds_ ) +
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

MultiArcTimeBiasParameter::MultiArcTimeBiasParameter( const EstimatebleParametersEnum parameterName,
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
    arcStartTimes_( arcStartTimes ), linkEndForTime_( linkEndForTime ), linkEnds_( linkEnds ), observableType_( observableType ),
    numberOfArcs_( static_cast< int >( arcStartTimes.size( ) ) )
{
    linkEndIndex_ =
            observation_models::getLinkEndIndicesForLinkEndTypeAtObservable( observableType_, linkEndForTime_, linkEnds_.size( ) ).at( 0 );
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
        throw std::runtime_error( "Error in " + getParameterTypeString( parameterName_.first ) + " of observable type " +
                                  observation_models::getObservableName( observableType_, linkEnds_.size( ) ) +
                                  " with link ends: " + observation_models::getLinkEndsString( linkEnds_ ) +
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

ConstantTimeDriftBiasParameter::ConstantTimeDriftBiasParameter( const EstimatebleParametersEnum parameterName,
                                                                const std::function< Eigen::VectorXd( ) > getCurrentBias,
                                                                const std::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
                                                                const int linkEndIndex,
                                                                const observation_models::LinkEnds linkEnds,
                                                                const observation_models::ObservableType observableType,
                                                                const double referenceEpoch,
                                                                const std::string& pointOnBodyId ):
    SingleArcObservationBiasParameter( parameterName, getCurrentBias, resetCurrentBias, linkEnds, observableType, pointOnBodyId ),
    linkEndIndex_( linkEndIndex ), referenceEpoch_( referenceEpoch )
{}

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
{}

std::vector< double > ArcWiseTimeDriftBiasParameter::getReferenceEpochs( )
{
    return referenceEpochs_;
}

}  // namespace estimatable_parameters

}  // namespace tudat
