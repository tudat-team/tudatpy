/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONBIASPARAMETER_H
#define TUDAT_OBSERVATIONBIASPARAMETER_H

#include <functional>
#include <iostream>
#include <vector>

#include "tudat/astro/observation_models/observationBias.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Utility function to create more precise secondary identification string of observation bias parameters.
std::string createObsBiasSecondaryIdentifier( const observation_models::ObservableType observableType,
                                              const observation_models::LinkEnds& linkEnds );

//! Common bias-function plumbing for all observation-bias estimatable parameters.
/*!
 * The wrapper owns the get/set callback pair while derived classes define
 * parameter sizing and metadata.
 */
template< typename BiasValueType >
class ObservationBiasFunctionWrapper : public EstimatableParameter< Eigen::VectorXd >
{
public:
    using GetBiasFunction = std::function< BiasValueType( ) >;

    using ResetBiasFunction = std::function< void( const BiasValueType& ) >;

    virtual ~ObservationBiasFunctionWrapper( ) {}

    //! Binds get/set callbacks during closure and applies any deferred value immediately.
    void setObservationBiasFunctions( const GetBiasFunction& getBiasFunction, const ResetBiasFunction& resetBiasFunction )
    {
        if( !( getBiasFunction_ == nullptr ) || !( resetBiasFunction_ == nullptr ) )
        {
            std::cerr << "Warning when resetting observation bias functions in estimation object, existing contents not empty" << std::endl;
        }

        getBiasFunction_ = getBiasFunction;
        resetBiasFunction_ = resetBiasFunction;

        // Apply a value that was set before closure, now that the reset callback is available.
        if( hasDeferredBiasValue_ && resetBiasFunction_ != nullptr )
        {
            resetBiasFunction_( deferredBiasValue_ );
            hasDeferredBiasValue_ = false;
        }
    }

protected:
    ObservationBiasFunctionWrapper( const EstimatebleParametersEnum parameterName,
                                    const std::string& associatedBody,
                                    const std::string& pointOnBodyId,
                                    const GetBiasFunction& getBiasFunction = GetBiasFunction( ),
                                    const ResetBiasFunction& resetBiasFunction = ResetBiasFunction( ) ):
        EstimatableParameter< Eigen::VectorXd >( parameterName, associatedBody, pointOnBodyId ), getBiasFunction_( getBiasFunction ),
        resetBiasFunction_( resetBiasFunction )
    {}

    bool biasFunctionsAreDefined( ) const
    {
        return ( getBiasFunction_ != nullptr && resetBiasFunction_ != nullptr );
    }

    //! Returns the bound getter callback for the linked observation-bias model.
    const GetBiasFunction& getBiasFunction( ) const
    {
        return getBiasFunction_;
    }

    //! Returns the bound reset callback for the linked observation-bias model.
    const ResetBiasFunction& resetBiasFunction( ) const
    {
        return resetBiasFunction_;
    }

    //! Stores a value that should be applied once closure creates reset callbacks.
    void deferBiasValue( const BiasValueType& biasValue )
    {
        deferredBiasValue_ = biasValue;
        hasDeferredBiasValue_ = true;
    }

    //! Indicates whether a pre-closure value assignment is waiting to be applied.
    bool hasDeferredBiasValue( ) const
    {
        return hasDeferredBiasValue_;
    }

    //! Returns the deferred value captured before closure was completed.
    const BiasValueType& getDeferredBiasValue( ) const
    {
        return deferredBiasValue_;
    }

    void clearDeferredBiasValue( )
    {
        hasDeferredBiasValue_ = false;
    }

    //! Writes to the linked model if available, otherwise stores the value for post-closure application.
    void resetOrDeferBiasValue( const BiasValueType& biasValue )
    {
        if( resetBiasFunction_ != nullptr )
        {
            resetBiasFunction_( biasValue );
            hasDeferredBiasValue_ = false;
        }
        else
        {
            deferBiasValue( biasValue );
        }
    }

private:
    //! Getter callback to read the current value from the underlying bias model.
    GetBiasFunction getBiasFunction_;

    //! Setter callback to write a new value to the underlying bias model.
    ResetBiasFunction resetBiasFunction_;

    //! Value captured by setParameterValue before closure links callbacks.
    BiasValueType deferredBiasValue_;

    //! True when deferredBiasValue_ should be applied at first successful closure.
    bool hasDeferredBiasValue_ = false;
};

//! Generic single-arc observation-bias parameter with Eigen::VectorXd callback interface.
class SingleArcObservationBiasParameter : public ObservationBiasFunctionWrapper< Eigen::VectorXd >
{
public:
    SingleArcObservationBiasParameter( const EstimatebleParametersEnum parameterName,
                                       const std::function< Eigen::VectorXd( ) > getCurrentBias,
                                       const std::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
                                       const observation_models::LinkEnds linkEnds,
                                       const observation_models::ObservableType observableType,
                                       const std::string& pointOnBodyId = "" );

    ~SingleArcObservationBiasParameter( ) {}

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( Eigen::VectorXd parameterValue );

    int getParameterSize( );

    void throwExceptionIfNotFullyDefined( );

    observation_models::LinkEnds getLinkEnds( );

    observation_models::ObservableType getObservableType( );

    std::string getParameterDescription( );

    //! Link geometry this estimatable bias parameter belongs to.
    observation_models::LinkEnds linkEnds_;

    //! Observable type this parameter applies to.
    observation_models::ObservableType observableType_;
};

//! Generic multi-arc observation-bias parameter with lookup-scheme support.
class MultiArcObservationBiasParameter : public ObservationBiasFunctionWrapper< std::vector< Eigen::VectorXd > >
{
public:
    MultiArcObservationBiasParameter( const EstimatebleParametersEnum parameterName,
                                      const std::vector< double > arcStartTimes,
                                      const std::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
                                      const std::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
                                      const int linkEndIndex,
                                      const observation_models::LinkEnds linkEnds,
                                      const observation_models::ObservableType observableType,
                                      const std::string& pointOnBodyId = "" );

    ~MultiArcObservationBiasParameter( ) {}

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( Eigen::VectorXd parameterValue );

    int getParameterSize( );

    void throwExceptionIfNotFullyDefined( );

    observation_models::LinkEnds getLinkEnds( );

    observation_models::ObservableType getObservableType( );

    std::string getParameterDescription( );

    std::vector< double > getArcStartTimes( );

    int getLinkEndIndex( );

    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( );

    void setLookupScheme( const std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme );

    //! Arc boundaries for the arc-wise bias values.
    const std::vector< double > arcStartTimes_;

    //! Link-end index used to select the current arc in lookupScheme_.
    int linkEndIndex_;

    //! Link geometry this estimatable bias parameter belongs to.
    observation_models::LinkEnds linkEnds_;

    //! Observable type this parameter applies to.
    observation_models::ObservableType observableType_;

    //! Dimension of a single-observable bias vector.
    int observableSize_;

    //! Number of arc-wise bias blocks.
    int numberOfArcs_;

    //! Runtime map from evaluation time to arc index.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

//! Non-parameter mixin with acceleration callback used by time-bias partials.
class TimeBiasParameterBase
{
public:
    TimeBiasParameterBase( ) = default;
    virtual ~TimeBiasParameterBase( ) {}

    void setBodyAccelerationFunction( const std::function< Eigen::VectorXd( const double ) > bodyAccelerationFunction );

    std::function< Eigen::VectorXd( const double ) > getBodyAccelerationFunction( );

protected:
    //! Body acceleration provider used in time-bias partial calculations.
    std::function< Eigen::VectorXd( const double ) > bodyAccelerationFunction_;
};

//! Generic single-arc time-bias parameter (constant time bias).
class SingleArcTimeBiasParameter : public ObservationBiasFunctionWrapper< Eigen::VectorXd >, public TimeBiasParameterBase
{
public:
    SingleArcTimeBiasParameter( const EstimatebleParametersEnum parameterName,
                                const std::function< Eigen::VectorXd( ) > getCurrentBias,
                                const std::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
                                const observation_models::LinkEndType linkEndForTime,
                                const observation_models::LinkEnds linkEnds,
                                const observation_models::ObservableType observableType,
                                const std::string& pointOnBodyId = "" );

    ~SingleArcTimeBiasParameter( ) {}

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( Eigen::VectorXd parameterValue );

    int getParameterSize( );

    void throwExceptionIfNotFullyDefined( );

    observation_models::LinkEnds getLinkEnds( );

    observation_models::LinkEndId getLinkEndId( );

    observation_models::LinkEndType getReferenceLinkEnd( );

    observation_models::ObservableType getObservableType( );

    std::string getParameterDescription( );

    int getLinkEndIndex( );

protected:
    //! Link-end type whose epoch defines the time-bias interpretation.
    observation_models::LinkEndType linkEndForTime_;

    //! Link geometry this estimatable time-bias parameter belongs to.
    observation_models::LinkEnds linkEnds_;

    //! Observable type this parameter applies to.
    observation_models::ObservableType observableType_;

    //! Cached index of linkEndForTime_ in the observable link-end ordering.
    int linkEndIndex_;
};

//! Generic multi-arc time-bias parameter (arc-wise time bias).
class MultiArcTimeBiasParameter : public ObservationBiasFunctionWrapper< std::vector< Eigen::VectorXd > >, public TimeBiasParameterBase
{
public:
    MultiArcTimeBiasParameter( const EstimatebleParametersEnum parameterName,
                               const std::vector< double > arcStartTimes,
                               const std::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
                               const std::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
                               const observation_models::LinkEndType linkEndForTime,
                               const observation_models::LinkEnds linkEnds,
                               const observation_models::ObservableType observableType,
                               const std::string& pointOnBodyId = "" );

    ~MultiArcTimeBiasParameter( ) {}

    Eigen::VectorXd getParameterValue( );

    void setParameterValue( Eigen::VectorXd parameterValue );

    int getParameterSize( );

    void throwExceptionIfNotFullyDefined( );

    observation_models::LinkEnds getLinkEnds( );

    observation_models::ObservableType getObservableType( );

    std::string getParameterDescription( );

    std::vector< double > getArcStartTimes( );

    int getLinkEndIndex( );

    std::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( );

    void setLookupScheme( const std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme );

    observation_models::LinkEndId getLinkEndId( );

    observation_models::LinkEndType getReferenceLinkEnd( );

protected:
    //! Arc boundaries for the arc-wise time-bias values.
    const std::vector< double > arcStartTimes_;

    //! Link-end type whose epoch defines the time-bias interpretation.
    observation_models::LinkEndType linkEndForTime_;

    //! Link geometry this estimatable time-bias parameter belongs to.
    observation_models::LinkEnds linkEnds_;

    //! Observable type this parameter applies to.
    observation_models::ObservableType observableType_;

    //! Number of arc-wise time-bias blocks.
    int numberOfArcs_;

    //! Cached index of linkEndForTime_ in the observable link-end ordering.
    int linkEndIndex_;

    //! Runtime map from evaluation time to arc index.
    std::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

//! Time-drift bias with specific reference-epoch metadata.
class ConstantTimeDriftBiasParameter : public SingleArcObservationBiasParameter
{
public:
    ConstantTimeDriftBiasParameter( const EstimatebleParametersEnum parameterName,
                                    const std::function< Eigen::VectorXd( ) > getCurrentBias,
                                    const std::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
                                    const int linkEndIndex,
                                    const observation_models::LinkEnds linkEnds,
                                    const observation_models::ObservableType observableType,
                                    const double referenceEpoch,
                                    const std::string& pointOnBodyId = "" );

    int getLinkEndIndex( );

    double getReferenceEpoch( );

private:
    //! Link-end index where the time-drift reference is applied.
    int linkEndIndex_;

    //! Reference epoch associated with the drift-bias parameterization.
    double referenceEpoch_;
};

//! Arc-wise time-drift bias with per-arc reference-epoch metadata.
class ArcWiseTimeDriftBiasParameter : public MultiArcObservationBiasParameter
{
public:
    ArcWiseTimeDriftBiasParameter( const EstimatebleParametersEnum parameterName,
                                   const std::vector< double > arcStartTimes,
                                   const std::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
                                   const std::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
                                   const int linkEndIndex,
                                   const observation_models::LinkEnds linkEnds,
                                   const observation_models::ObservableType observableType,
                                   const std::vector< double > referenceEpochs,
                                   const std::string& pointOnBodyId = "" );

    std::vector< double > getReferenceEpochs( );

private:
    //! Per-arc reference epochs matching arcStartTimes_.
    std::vector< double > referenceEpochs_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONBIASPARAMETER_H
