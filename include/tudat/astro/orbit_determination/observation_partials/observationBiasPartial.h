/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONBIASPARTIAL_H
#define TUDAT_OBSERVATIONBIASPARTIAL_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/observation_partials/directObservationPartial.h"

namespace tudat
{

namespace observation_partials
{

template< int ObservationSize >
class TimeBiasPartial
{
public:

    TimeBiasPartial(
        const observation_models::LinkEnds linkEnds,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_partials::DirectObservationPartial< ObservationSize > > transmitterWrtLinkEndStatePartial,
        const std::shared_ptr< observation_partials::DirectObservationPartial< ObservationSize > > receiverWrtLinkEndStatePartial,
        const std::function< Eigen::Vector3d( const double ) > transmitterAccelerationFunction,
        const std::function< Eigen::Vector3d( const double ) > receiverAccelerationFunction ):
            linkEnds_( linkEnds ), observableType_( observableType ), referenceLinkEnd_( referenceLinkEnd ),
            transmitterWrtLinkEndStatePartial_( transmitterWrtLinkEndStatePartial ),
            receiverWrtLinkEndStatePartial_( receiverWrtLinkEndStatePartial ),
            transmitterAccelerationFunction_( transmitterAccelerationFunction ),
            receiverAccelerationFunction_( receiverAccelerationFunction ){ }

    std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > getObservationPartialWrtObservationTime(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
        const Eigen::Matrix< double, ObservationSize, 1 >& currentObservation =
        Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) )
    {
        if( referenceLinkEnd_ != linkEndOfFixedTime )
        {
            return std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double >( );
        }
        else
        {
            double referenceTime = TUDAT_NAN;
            if( linkEndOfFixedTime == observation_models::transmitter )
            {
                transmitterWrtLinkEndStatePartial_->setFixLinkEndTime( true );
                receiverWrtLinkEndStatePartial_->setFixLinkEndTime( false );
                referenceTime = times.at( 0 );
            }
            else
            {
                transmitterWrtLinkEndStatePartial_->setFixLinkEndTime( false );
                receiverWrtLinkEndStatePartial_->setFixLinkEndTime( true );
                referenceTime = times.at( 1 );
            }

            std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > transmitterStatePartials =
                transmitterWrtLinkEndStatePartial_->calculatePartial(
                    states, times, linkEndOfFixedTime, ancillarySettings, currentObservation );
            std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > receiverStatePartials =
                receiverWrtLinkEndStatePartial_->calculatePartial(
                    states, times, linkEndOfFixedTime, ancillarySettings, currentObservation );

            transmitterWrtLinkEndStatePartial_->setFixLinkEndTime( false );
            receiverWrtLinkEndStatePartial_->setFixLinkEndTime( false );

            Eigen::Vector6d transmitterStateDerivative = Eigen::Vector6d::Zero( );
            transmitterStateDerivative.segment( 0, 3 ) = states.at( 0 ).segment( 3, 3 );
            transmitterStateDerivative.segment( 3, 3 ) = transmitterAccelerationFunction_( times.at( 0 ) );

            Eigen::Vector6d receiverStateDerivative = Eigen::Vector6d::Zero( );
            receiverStateDerivative.segment( 0, 3 ) = states.at( 1 ).segment( 3, 3 );
            receiverStateDerivative.segment( 3, 3 ) = receiverAccelerationFunction_( times.at( 1 ) );

            return std::make_pair( -( transmitterStatePartials.at( 0 ).first * transmitterStateDerivative +
                                      receiverStatePartials.at( 0 ).first * receiverStateDerivative ),
                                   referenceTime );
        }
    }

    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }


protected:

    observation_models::LinkEnds linkEnds_;

    observation_models::ObservableType observableType_;

    observation_models::LinkEndType referenceLinkEnd_;

    std::shared_ptr< observation_partials::DirectObservationPartial< ObservationSize > > transmitterWrtLinkEndStatePartial_;

    std::shared_ptr< observation_partials::DirectObservationPartial< ObservationSize > > receiverWrtLinkEndStatePartial_;

    std::function< Eigen::Vector3d( const double ) > transmitterAccelerationFunction_;

    std::function< Eigen::Vector3d( const double ) > receiverAccelerationFunction_;
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
    ObservationPartialWrtConstantTimeBias( const std::shared_ptr< TimeBiasPartial< ObservationSize > > timeBiasPartial ):
            ObservationPartial< ObservationSize >(
                    std::make_pair( estimatable_parameters::constant_time_observation_bias, timeBiasPartial->getLinkEnds( ).begin( )->second.getDualStringLinkEnd( ) ) ),
            timeBiasPartial_( timeBiasPartial )
    {
        if( timeBiasPartial == nullptr )
        {
            throw std::runtime_error( "Error, time bias partial not yet created." );
        }
    }

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
        return { timeBiasPartial_->getObservationPartialWrtObservationTime(
            states, times, linkEndOfFixedTime, ancillarySettings, currentObservation ) };
    }

private:

    std::shared_ptr< TimeBiasPartial< ObservationSize > > timeBiasPartial_;
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
    ObservationPartialWrtArcWiseTimeBias( const std::shared_ptr< TimeBiasPartial< ObservationSize > > timeBiasPartial,
                                          const std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme,
                                          const int linkEndIndex,
                                          const int numberOfArcs ):
            ObservationPartial< ObservationSize >( std::make_pair( estimatable_parameters::arc_wise_time_observation_bias,
                                                                   timeBiasPartial->getLinkEnds( ).begin( )->second.getDualStringLinkEnd( ) ) ),
            timeBiasPartial_( timeBiasPartial ), arcLookupScheme_( arcLookupScheme ),
            linkEndIndex_( linkEndIndex ), numberOfArcs_( numberOfArcs )
    {
        totalPartial_ = Eigen::MatrixXd::Zero( ObservationSize, numberOfArcs_ );
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
        totalPartial_.setZero( );
        int currentIndex = arcLookupScheme_->findNearestLowerNeighbour( times.at( 1 ) );

        std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > timeBiasPartialPair = timeBiasPartial_->getObservationPartialWrtObservationTime(
            states, times, linkEndOfFixedTime, ancillarySettings, currentObservation );
        totalPartial_.block( 0, currentIndex, ObservationSize, 1 ) = timeBiasPartialPair.first;
        return { std::make_pair( totalPartial_, timeBiasPartialPair.second ) };
    }

private:

    std::shared_ptr< TimeBiasPartial< ObservationSize > > timeBiasPartial_;

    //! Object used to determine the index from observationBiases_ to be used, based on the current time.
    std::shared_ptr< interpolators::LookUpScheme< double > > arcLookupScheme_;

    //! Link end index from which the 'current time' is determined
    int linkEndIndex_;

    //! Number of arcs for which biases are defined
    int numberOfArcs_;

    //! Pre-allocated partial vector
    Eigen::Matrix< double, ObservationSize, Eigen::Dynamic > totalPartial_;

};




}

}

#endif // OBSERVATIONPARTIAL_H
