/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONPARTIAL_H
#define TUDAT_OBSERVATIONPARTIAL_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/observationBiasParameter.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace observation_partials
{

//! Base class for scaling position partials to observable partials.
/*!
 *  Base class for scaling position partials to observable partials. For observables computed from the three-dimensional
 *  positions of reference points, the partial derivatives are computed from the partials of these positions. Each
 *  observation model requires a specific derived class of PositionPartialScaling, that computes how position partials
 *  are converted to ('scaled to') observable partials
 */
class PositionPartialScaling
{
public:

    //! Destructor
    virtual ~PositionPartialScaling( ){ }

    //! Update the scaling object to the current times and states (pure virtual).
    /*!
     *  Update the scaling object to the current times and states (pure virtual).
     *  \param linkEndStates List of states at each link end during observation.
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     */
    virtual void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                         const std::vector< double >& times,
                         const observation_models::LinkEndType fixedLinkEnd,
                         const Eigen::VectorXd currentObservation ) = 0;


};

template< int ObservationSize >
class DirectPositionPartialScaling: public PositionPartialScaling
{
public:

    DirectPositionPartialScaling( const observation_models::ObservableType observableType ):
        observableType_( observableType )
    { }

    virtual ~DirectPositionPartialScaling( ){ }

    virtual observation_models::LinkEndType getCurrentLinkEndType( ) = 0;

    virtual Eigen::Matrix< double, ObservationSize, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType ) = 0;

    virtual Eigen::Matrix< double, ObservationSize, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return Eigen::MatrixXd::Zero( ObservationSize, 3 );
    }

    virtual Eigen::Matrix< double, ObservationSize, 1 > getLightTimePartialScalingFactor( ) = 0;

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    virtual bool isVelocityScalingNonZero( )
    {
        return false;
    }

    virtual std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > getLinkIndependentPartials(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterType )
    {
        if( useLinkIndependentPartials( ) == false )
        {
            throw std::runtime_error( "Error when getting link independent partials of observable type " +
                                      observation_models::getObservableName( observableType_) +
                                      ", no such model is defined." );
        }
        else
        {
            throw std::runtime_error( "Error when getting link independent partials of observable type " +
                                      observation_models::getObservableName( observableType_) +
                                      ", derived class model is not implemented." );
        }
    }

    virtual bool useLinkIndependentPartials( )
    {
        return false;
    }

protected:

    observation_models::ObservableType observableType_;

    bool doesVelocityScalingFactorExist_;
};

//! Base class of partial derivative of an observable w.r.t. an estimated parameter.
/*!
 *  Base class of partial derivative of an observable w.r.t. an estimated parameter.
 *  Calculated partial function must be defined in each derived class for calculation of specific derivative.
 */
template< int ObservationSize >
class ObservationPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param parameterIdentifier Parameter id of for specifc parameter of which partial is computed by object.
     */
    ObservationPartial( const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        parameterIdentifier_( parameterIdentifier ){ }

    //! Virtual destructor.
    /*!
     *  Virtual destrcutor.
     */
    virtual ~ObservationPartial( ) { }

    //! Pure virtual function to calculate the observation partial(s) at required time(s) and state(s)
    /*!
     *  Pure virtual function to calculate the observation partial(s) at required time(s) and state(s). States and times
     *  are typically obtained from evaluation of associated observation model.
     *  Derived class functions implement this for a specific observable.
     *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
     *  \param times Link end times.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes).
     *  \return Vector of pairs containing partial values and associated times.
     */
    virtual std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) ) = 0;

    //! Function to get parameter id of for specifc parameter of which partial is computed by object.
    /*!
     * Function to get parameter id of for specifc parameter of which partial is computed by object.
     * \return Parameter id of for specifc parameter of which partial is computed by object.
     */
    estimatable_parameters::EstimatebleParameterIdentifier getParameterIdentifier( )
    {
        return parameterIdentifier_;
    }


protected:

    //! Parameter id of for specifc parameter of which partial is computed by object.
    estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier_;
};


//! Class for computing the derivative of any observable w.r.t. a constant absolute observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a constant absolute observation bias. Note that this partial is
 *  distinct from most other ObservationPartial partial derived classes, as its implementation is based on the parameter
 *  (constant observation bias), not the type of observable: the implementation is identical for each observable.
 */
template< int ObservationSize >
class ObservationPartialWrtConstantAbsoluteBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     */
    ObservationPartialWrtConstantAbsoluteBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::constant_additive_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds )
    {
        // Compute partial (vector of ObservationSize with 1.0 entries).
        constantPartial_ = Eigen::Matrix< double, ObservationSize, ObservationSize >::Identity( );
    }

    //! Destructor
    ~ObservationPartialWrtConstantAbsoluteBias( ){ }

    //! Function to calculate the observation partial w.r.t. constant absolute bias
    /*!
     *  Function to calculate the observation partial w.r.t. constant absolute bias. Note that output is independent of input.
     *  Associated time defined at times[ 0 ].
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
        return { std::make_pair( constantPartial_, times.at( 0 ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Observation partial: constant for all conditions.
    Eigen::Matrix< double, ObservationSize, ObservationSize > constantPartial_;

};


//! Class for computing the derivative of any observable w.r.t. an arc-wise constant absolute observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a n arc-wiseconstant absolute observation bias. Note that this
 *  partial is distinct from most other ObservationPartial partial derived classes, as its implementation is based on the
 *  parameter (arc-wise constant observation bias), not the type of observable: the implementation is identical for each
 *  observable.
 */
template< int ObservationSize >
class ObservationPartialWrtArcWiseAbsoluteBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param arcLookupScheme Object used to determine the index from observationBiases_ to be used, based on the current time.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     * \param numberOfArcs Number of arcs for which biases are defined
     */
    ObservationPartialWrtArcWiseAbsoluteBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds,
                                               const std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme,
                                               const int linkEndIndex,
                                               const int numberOfArcs ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::arcwise_constant_additive_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds ), arcLookupScheme_( arcLookupScheme ),
        linkEndIndex_( linkEndIndex ), numberOfArcs_( numberOfArcs )
    {
        constantPartial_ = Eigen::Matrix< double, ObservationSize, ObservationSize >::Identity( );
        totalPartial_ = Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >::Zero(
                    ObservationSize, ObservationSize * numberOfArcs_ );
    }

    //! Destructor
    ~ObservationPartialWrtArcWiseAbsoluteBias( ){ }

    //! Function to calculate the observation partial w.r.t. arc-wise constant absolute bias
    /*!
     *  Function to calculate the observation partial w.r.t. arc-wise constant absolute bias. Note that output is independent of
     *  input. Associated time defined by linkEndIndex_.
     *  \param states Link end states (unused).
     *  \param times Link end times.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
//        std::cout << "TEST CALCULATE PARTIALS W.R.T. ARC WISE CONSTANT BIAS" << "\n\n";
        int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );


        totalPartial_.setZero( );
        if( arcLookupScheme_->getMinimumValue( ) <= times.at( linkEndIndex_ ) )
        {
            int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );
            totalPartial_.block( 0, currentIndex * ObservationSize, ObservationSize, ObservationSize ) = constantPartial_;
        }

        return { std::make_pair( totalPartial_, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Number of arcs for which biases are defined
    int numberOfArcs_;

       //! Observation partial for singe arc: constant for all conditions.
    Eigen::Matrix< double, ObservationSize, ObservationSize > constantPartial_;

    //! Pre-allocated partial vector
    Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > totalPartial_;

};

//! Class for computing the derivative of any observable w.r.t. a constant relative observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a constant relative observation bias. Note that this partial is
 *  distinct from most other ObservationPartial partial derived classes, as its implementation is based on the parameter
 *  (constant observation bias), not the type of observable: the implementation is identical for each observable.
 */
template< int ObservationSize >
class ObservationPartialWrtConstantRelativeBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     */
    ObservationPartialWrtConstantRelativeBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::constant_additive_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds )
    {  }

    //! Destructor
    ~ObservationPartialWrtConstantRelativeBias( ){ }

    //! Function to calculate the observation partial w.r.t. constant relative bias
    /*!
     *  Function to calculate the observation partial w.r.t. constant bias. Note that output is independent of input times and
     *  states. Associated time defined at times[ 0 ].
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed.
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        return { std::make_pair( currentObservation, times.at( 0 ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;
};

//! Class for computing the derivative of any observable w.r.t. an arc-wise constant relative observation bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a n arc-wiseconstant relative observation bias. Note that this
 *  partial is distinct from most other ObservationPartial partial derived classes, as its implementation is based on the
 *  parameter (arc-wise constant observation bias), not the type of observable: the implementation is identical for each
 *  observable.
 */
template< int ObservationSize >
class ObservationPartialWrtArcWiseRelativeBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param arcLookupScheme Object used to determine the index from observationBiases_ to be used, based on the current time.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     * \param numberOfArcs Number of arcs for which biases are defined
     */
    ObservationPartialWrtArcWiseRelativeBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds,
                                               const std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme,
                                               const int linkEndIndex,
                                               const int numberOfArcs ):
        ObservationPartial< ObservationSize >(
            std::make_pair( estimatable_parameters::arcwise_constant_relative_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( ) ) ),
        observableType_( observableType ), linkEnds_( linkEnds ), arcLookupScheme_( arcLookupScheme ),
        linkEndIndex_( linkEndIndex ), numberOfArcs_( numberOfArcs )
    {
       totalPartial_ = Eigen::VectorXd::Zero( ObservationSize * numberOfArcs_ );
    }

    //! Destructor
    ~ObservationPartialWrtArcWiseRelativeBias( ){ }

    //! Function to calculate the observation partial w.r.t. arc-wise constant relative bias
    /*!
     *  Function to calculate the observation partial w.r.t. arc-wise constant relative bias. Note that output is independent of
     *  input. Associated time defined by linkEndIndex_.
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
        totalPartial_.setZero( );

        if( arcLookupScheme_->getMinimumValue( ) <= times.at( linkEndIndex_ ) )
        {
            int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );
            totalPartial_.segment( currentIndex * ObservationSize, ObservationSize ) = currentObservation;
        }
        return { std::make_pair( totalPartial_, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Number of arcs for which biases are defined
    int numberOfArcs_;

    //! Pre-allocated partial vector
    Eigen::VectorXd totalPartial_;

};

//! Class for computing the derivative of any observable w.r.t. a constant time drift bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a constant time drift bias. Note that this partial is
 *  distinct from most other ObservationPartial partial derived classes, as its implementation is based on the parameter
 *  (constant time drift bias), not the type of observable: the implementation is identical for each observable.
 */
template< int ObservationSize >
class ObservationPartialWrtConstantTimeDriftBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     * \param referenceEpoch Reference epoch at which the time drift is initialised.
     */
    ObservationPartialWrtConstantTimeDriftBias( const observation_models::ObservableType observableType,
                                                const observation_models::LinkEnds& linkEnds,
                                                const int linkEndIndex,
                                                const double referenceEpoch ):
            ObservationPartial< ObservationSize >(
                    std::make_pair( estimatable_parameters::constant_time_drift_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( )  ) ),
            observableType_( observableType ), linkEnds_( linkEnds ), linkEndIndex_( linkEndIndex ), referenceEpoch_( referenceEpoch )
    {  }

    //! Destructor
    ~ObservationPartialWrtConstantTimeDriftBias( ){ }

    //! Function to calculate the observation partial w.r.t. time drift bias
    /*!
     *  Function to calculate the observation partial w.r.t. time drift bias. Note that output is dependent on input times but independent of
     *  states. Associated time defined at times[ 0 ].
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed.
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        Eigen::Matrix< double, ObservationSize, 1 > observationTime;
        for ( unsigned int i = 0 ; i < ObservationSize ; i++ )
        {
            observationTime( i, 0 ) = times.at( linkEndIndex_ ) - referenceEpoch_;
        }
        return { std::make_pair( observationTime, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Reference epoch at which the time drift is initialised.
    double referenceEpoch_;
};


//! Class for computing the derivative of any observable w.r.t. an arc-wise time drift bias
/*!
*  Class for computing the derivative of any observable w.r.t. an arc-wise time drift bias. Note that this
*  partial is distinct from most other ObservationPartial partial derived classes, as its implementation is based on the
*  parameter (arc-wise time drift bias), not the type of observable: the implementation is identical for each
*  observable.
*/
template< int ObservationSize >
class ObservationPartialWrtArcWiseTimeDriftBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param arcLookupScheme Object used to determine the index from observationBiases_ to be used, based on the current time.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     * \param numberOfArcs Number of arcs for which biases are defined
     * \param referenceEpochs Reference epochs (for each arc) at which the time drifts are initialised
     */
    ObservationPartialWrtArcWiseTimeDriftBias( const observation_models::ObservableType observableType,
                                          const observation_models::LinkEnds& linkEnds,
                                          const std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme,
                                          const int linkEndIndex,
                                          const int numberOfArcs,
                                          const std::vector< double > referenceEpochs ):
            ObservationPartial< ObservationSize >( std::make_pair( estimatable_parameters::arc_wise_time_drift_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( )  ) ),
            observableType_( observableType ), linkEnds_( linkEnds ), arcLookupScheme_( arcLookupScheme ),
            linkEndIndex_( linkEndIndex ), numberOfArcs_( numberOfArcs ), referenceEpochs_( referenceEpochs )
    {
        totalPartial_ = Eigen::VectorXd::Zero( ObservationSize * numberOfArcs_ );
    }

    //! Destructor
    ~ObservationPartialWrtArcWiseTimeDriftBias( ){ }

    //! Function to calculate the observation partial w.r.t. arc-wise time drift bias
    /*!
     *  Function to calculate the observation partial w.r.t. arc-wise time drift bias. Note that output is independent of
     *  input. Associated time defined by linkEndIndex_.
     *  \param states Link end states (unused).
     *  \param times Link end times.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
        int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );

        Eigen::Matrix< double, ObservationSize, 1 > observationTime;
        for ( unsigned int i = 0 ; i < ObservationSize ; i++ )
        {
            observationTime( i, 0 ) = times.at( linkEndIndex_ ) - referenceEpochs_.at( currentIndex );
        }

        totalPartial_.setZero( );
        totalPartial_.segment( currentIndex * ObservationSize, ObservationSize ) = observationTime;

        return { std::make_pair( totalPartial_, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Number of arcs for which biases are defined
    int numberOfArcs_;

    //! Pre-allocated partial vector
    Eigen::VectorXd totalPartial_;

    //! Reference epochs (per arc) at which the time drift biases are initialised.
    std::vector< double > referenceEpochs_;

};

//! Class for computing the derivative of any observable w.r.t. a constant time bias
/*!
 *  Class for computing the derivative of any observable w.r.t. a constant time bias. Note that this partial is
 *  distinct from most other ObservationPartial partial derived classes, as its implementation is based on the parameter
 *  (constant time bias), not the type of observable: the implementation is identical for each observable.
 */
template< int ObservationSize >
class ObservationPartialWrtConstantTimeBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     */
    ObservationPartialWrtConstantTimeBias( const observation_models::ObservableType observableType,
                                           const observation_models::LinkEnds& linkEnds,
                                           const int linkEndIndex ):
            ObservationPartial< ObservationSize >(
                    std::make_pair( estimatable_parameters::constant_time_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( )  ) ),
            observableType_( observableType ), linkEnds_( linkEnds ), linkEndIndex_( linkEndIndex )
    {  }

    //! Destructor
    ~ObservationPartialWrtConstantTimeBias( ){ }

    //! Function to calculate the observation partial w.r.t. time bias
    /*!
     *  Function to calculate the observation partial w.r.t. time bias. Note that output is dependent on input times but independent of
     *  states. Associated time defined at times[ 0 ].
     *  \param states Link end states (unused).
     *  \param times Link end times  (unused).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable  (unused).
     *  \param currentObservation Value of the observation for which the partial is to be computed.
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        Eigen::Matrix< double, ObservationSize, 1 > observationTime;
        for ( unsigned int i = 0 ; i < ObservationSize ; i++ )
        {
            observationTime( i, 0 ) = - 1.0;
        }
        return { std::make_pair( observationTime, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;
};


//! Class for computing the derivative of any observable w.r.t. an arc-wise time bias
/*!
*  Class for computing the derivative of any observable w.r.t. an arc-wise time bias. Note that this
*  partial is distinct from most other ObservationPartial partial derived classes, as its implementation is based on the
*  parameter (arc-wise time bias), not the type of observable: the implementation is identical for each
*  observable.
*/
template< int ObservationSize >
class ObservationPartialWrtArcWiseTimeBias: public ObservationPartial< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type for which the bias is active.
     * \param linkEnds Observation link ends for which the bias is active.
     * \param arcLookupScheme Object used to determine the index from observationBiases_ to be used, based on the current time.
     * \param linkEndIndex Link end index from which the 'current time' is determined
     * \param numberOfArcs Number of arcs for which biases are defined
     */
    ObservationPartialWrtArcWiseTimeBias( const observation_models::ObservableType observableType,
                                               const observation_models::LinkEnds& linkEnds,
                                               const std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme,
                                               const int linkEndIndex,
                                               const int numberOfArcs ):
            ObservationPartial< ObservationSize >( std::make_pair( estimatable_parameters::arc_wise_time_observation_bias, linkEnds.begin( )->second.getDualStringLinkEnd( )  ) ),
            observableType_( observableType ), linkEnds_( linkEnds ), arcLookupScheme_( arcLookupScheme ),
            linkEndIndex_( linkEndIndex ), numberOfArcs_( numberOfArcs )
    {
        totalPartial_ = Eigen::VectorXd::Zero( 1 * numberOfArcs_ );
    }

    //! Destructor
    ~ObservationPartialWrtArcWiseTimeBias( ){ }

    //! Function to calculate the observation partial w.r.t. arc-wise time bias
    /*!
     *  Function to calculate the observation partial w.r.t. arc-wise time bias. Note that output is independent of
     *  input. Associated time defined by linkEndIndex_.
     *  \param states Link end states (unused).
     *  \param times Link end times.
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed  (unused).
     *  \return Vector of pairs containing partial values and associated times.
     */
    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
                Eigen::Matrix< double, ObservationSize, 1 >::Zero( ) )
    {
        int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( linkEndIndex_ ) );

        Eigen::Matrix< double, ObservationSize, 1 > observationTime;
        for ( unsigned int i = 0 ; i < ObservationSize ; i++ )
        {
            observationTime( i, 0 ) = - 1.0;
        }

        totalPartial_.setZero( );
        totalPartial_.segment( currentIndex * ObservationSize, ObservationSize ) = observationTime;

        return { std::make_pair( totalPartial_, times.at( linkEndIndex_ ) ) };
    }

private:

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    //!  Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Number of arcs for which biases are defined
    int numberOfArcs_;

    //! Pre-allocated partial vector
    Eigen::VectorXd totalPartial_;

};


//extern template class ObservationPartial< 1 >;
//extern template class ObservationPartial< 2 >;
//extern template class ObservationPartial< 3 >;

//extern template class ObservationPartialWrtConstantAbsoluteBias< 1 >;
//extern template class ObservationPartialWrtConstantAbsoluteBias< 2 >;
//extern template class ObservationPartialWrtConstantAbsoluteBias< 3 >;

//extern template class ObservationPartialWrtArcWiseAbsoluteBias< 1 >;
//extern template class ObservationPartialWrtArcWiseAbsoluteBias< 2 >;
//extern template class ObservationPartialWrtArcWiseAbsoluteBias< 3 >;

//extern template class ObservationPartialWrtConstantRelativeBias< 1 >;
//extern template class ObservationPartialWrtConstantRelativeBias< 2 >;
//extern template class ObservationPartialWrtConstantRelativeBias< 3 >;

//extern template class ObservationPartialWrtArcWiseRelativeBias< 1 >;
//extern template class ObservationPartialWrtArcWiseRelativeBias< 2 >;
//extern template class ObservationPartialWrtArcWiseRelativeBias< 3 >;




//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 1.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > > SingleLinkObservationPartialList;

//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 2.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 2 > > > SingleLinkObservationTwoPartialList;

//! Typedef for map of observation partials.
/*!
 *  Typedef for map of observation partials, for an observable of size 3.
 *  The key of the map is a pair denoting the start index of the parameter
 *  in the list of estimated parameters and the number of indices of the parameter
 *  wrt which the current partial (corresponding value in map) is taken.
 */
typedef std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 3 > > > SingleLinkObservationThreePartialList;


//! Function to create partials of observation w.r.t. a link property.
/*!
 *  Function to create partials of observation w.r.t. a link property, e.g. a parameter that does not influence either link end's
 *  dynamics, only the observable itself, such as observation biases and clock parameters.
 *  \param linkEnds Link ends of observable for which partial is to be made.
 *  \param observableType Type of observable for which partial is to be made.
 *  \param parameterToEstimate Parameter w.r.t. which the partial is to be taken
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Object that computes the partial of the observation w.r.t. parameterToEstimate (nullptr if no dependency).
 */
template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtLinkProperty(
        const observation_models::LinkEnds& linkEnds,
        const observation_models::ObservableType observableType,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterToEstimate,
        const bool useBiasPartials = true )
{
//    std::cout << "observableType: " << observableType << "\n\n";

    std::shared_ptr< ObservationPartial< ObservationSize > > observationPartial;

    // Check parameter type
    switch( parameterToEstimate->getParameterName( ).first )
    {
    case estimatable_parameters::constant_additive_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > constantBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                        parameterToEstimate );
            if( constantBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == constantBias->getLinkEnds( ) && observableType == constantBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtConstantAbsoluteBias< ObservationSize > >(
                                observableType, linkEnds );
                }
            }
        }
        break;
    }
    case estimatable_parameters::arcwise_constant_additive_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseObservationBiasParameter > arcwiseBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseObservationBiasParameter >(
                        parameterToEstimate );
            if( arcwiseBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. arcwise observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == arcwiseBias->getLinkEnds( ) && observableType == arcwiseBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseAbsoluteBias< ObservationSize > >(
                                observableType, linkEnds,
                                arcwiseBias->getLookupScheme( ),
                                arcwiseBias->getLinkEndIndex( ),
                                arcwiseBias->getArcStartTimes( ).size( ) );
                }
            }
        }
        break;
    }
    case estimatable_parameters::constant_relative_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ConstantObservationBiasParameter > constantBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantObservationBiasParameter >(
                        parameterToEstimate );
            if( constantBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == constantBias->getLinkEnds( ) && observableType == constantBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtConstantRelativeBias< ObservationSize > >(
                                observableType, linkEnds );

                }
            }
        }
        break;
    }
    case estimatable_parameters::arcwise_constant_relative_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseObservationBiasParameter > arcwiseBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseObservationBiasParameter >(
                        parameterToEstimate );
            if( arcwiseBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. arcwise relative observation bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == arcwiseBias->getLinkEnds( ) && observableType == arcwiseBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseRelativeBias< ObservationSize > >(
                                observableType, linkEnds,
                                arcwiseBias->getLookupScheme( ),
                                arcwiseBias->getLinkEndIndex( ),
                                arcwiseBias->getArcStartTimes( ).size( ) );
                }
            }
        }
        break;
    }
    case estimatable_parameters::constant_time_drift_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ConstantTimeDriftBiasParameter > constantTimeDriftBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantTimeDriftBiasParameter >(
                            parameterToEstimate );
            if( constantTimeDriftBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. time drift bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == constantTimeDriftBias->getLinkEnds( ) && observableType == constantTimeDriftBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtConstantTimeDriftBias< ObservationSize > >(
                            observableType, linkEnds, constantTimeDriftBias->getLinkEndIndex( ), constantTimeDriftBias->getReferenceEpoch( ) );
                }
            }
        }
        break;
    }
    case estimatable_parameters::arc_wise_time_drift_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseTimeDriftBiasParameter > arcwiseTimeDriftBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseTimeDriftBiasParameter >( parameterToEstimate );
            if ( arcwiseTimeDriftBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. arcwise time drift bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if ( linkEnds == arcwiseTimeDriftBias->getLinkEnds( ) && observableType == arcwiseTimeDriftBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseTimeDriftBias< ObservationSize > >(
                            observableType, linkEnds,
                            arcwiseTimeDriftBias->getLookupScheme( ),
                            arcwiseTimeDriftBias->getLinkEndIndex( ),
                            arcwiseTimeDriftBias->getArcStartTimes( ).size( ),
                            arcwiseTimeDriftBias->getReferenceEpochs( ) );
                }
            }
        }
        break;
    }
    case estimatable_parameters::constant_time_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ConstantTimeBiasParameter > constantTimeBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ConstantTimeBiasParameter >( parameterToEstimate );
            if( constantTimeBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. time bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if( linkEnds == constantTimeBias->getLinkEnds( ) && observableType == constantTimeBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtConstantTimeBias< ObservationSize > >(
                            observableType, linkEnds, constantTimeBias->getLinkEndIndex( ) );
                }
            }
        }
        break;
    }
    case estimatable_parameters::arc_wise_time_observation_bias:
    {
        if( useBiasPartials )
        {
            // Check input consistency
            std::shared_ptr< estimatable_parameters::ArcWiseTimeBiasParameter > arcwiseTimeBias =
                    std::dynamic_pointer_cast< estimatable_parameters::ArcWiseTimeBiasParameter >( parameterToEstimate );
            if ( arcwiseTimeBias == nullptr )
            {
                throw std::runtime_error( "Error when making partial w.r.t. arcwise time bias, type is inconsistent" );
            }
            else
            {
                // Check dependency between parameter and link properties.
                if ( linkEnds == arcwiseTimeBias->getLinkEnds( ) && observableType == arcwiseTimeBias->getObservableType( ) )
                {
                    observationPartial = std::make_shared< ObservationPartialWrtArcWiseTimeBias< ObservationSize > >(
                            observableType, linkEnds,
                            arcwiseTimeBias->getLookupScheme( ),
                            arcwiseTimeBias->getLinkEndIndex( ),
                            arcwiseTimeBias->getArcStartTimes( ).size( ) );
                }
            }
        }
        break;
    }
    default:
        break;
    }

    return observationPartial;
}

}

}

#endif // OBSERVATIONPARTIAL_H
