/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NWAYRANGEOBSERVATIONMODEL_H
#define TUDAT_NWAYRANGEOBSERVATIONMODEL_H

#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
namespace tudat
{

namespace observation_models
{

//! Class for simulating n-way range observables.
/*!
 *  Class for simulating n-way range observations. It connects an arbitrary number of links by rane observables to obtain the
 *  n-way range. In most practical cases (e.g. DSN radio ranging, SLR), n will be equal to 2. Note that here the number of 'ways'
 *  represents the number of legs that the signal travels along. As a result, for both DSN 2-way and 3-way observables, n equals
 *  2 in this model. The difference is that the first and last link end will be the same for the former case, and different for
 *  the latter case. The retransmission of a signal at the intermediate link ends can be done with a (negative or positive) delay.
 */
template< typename ObservationScalarType = double,
          typename TimeType = double >
class NWayRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:    
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculators List of objects to compute the light-times (including any corrections w.r.t. Euclidean case)
     *  for each leg of the n-way range. First entry starts at transmitter; last entry is to receiver.
     *  \param retransmissionDelays Function that returns the list of retransmission delays as a function of observation time.
     *  The retransmission delays represent the time difference between the reception of a singal by an intermediate link end,
     *  and the retransmission to the subsequent link end. By default, this function is empty, in which case no retransmission
     *  delays are used.
     *  \param observationBiasCalculator Object for calculating (system-dependent) errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    NWayRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::vector< std::shared_ptr< observation_models::LightTimeCalculator
            < ObservationScalarType, TimeType > > > lightTimeCalculators,
            const std::function< std::vector< double >( const double ) > retransmissionDelays =
            std::function< std::vector< double >( const double ) >( ),
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_range, linkEnds, observationBiasCalculator ),
        lightTimeCalculators_( lightTimeCalculators ), retransmissionDelays_( retransmissionDelays )
    {
        numberOfLinks_ = lightTimeCalculators_.size( );
        numberOfLinkEnds_ = numberOfLinks_ + 1;
    }

    //! Destructor
    ~NWayRangeObservationModel( ){ }

    //! Function to compute n-way range observable without any corrections.
    /*!
     *  Function to compute n-way range  observable without any corrections, i.e. the true physical range as computed
     *  from the defined link ends. The time argument can be at any of the link ends
     *  involved in the onbservation (including the intermediate link ends) by the linkEndAssociatedWithTime input).
     *  In the case where the reference link end is an intermediate link end, the inpit time denotes the reception time
     *  of the signal at this station (which need not be the same as its retransmission time).
     *  Note that this observable does include light-time corrections, which represent physically true corrections. It does not
     *  include e.g. system-dependent measurement errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal n-way range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        // Initialize total light-time/single-leg light-time
        ObservationScalarType totalLightTime =
                mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );
        ObservationScalarType currentLightTime;
        StateType currentReceiverStateOutput, currentTransmitterStateOutput;

        // Resize link-end states/times
        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndTimes.resize( 2 * ( numberOfLinkEnds_ - 1 ) );
        linkEndStates.resize( 2 * ( numberOfLinkEnds_ - 1 ) );

        // Retrieve retransmission delays
        if( !( retransmissionDelays_ == nullptr ) )
        {
            currentRetransmissionDelays_ = retransmissionDelays_( time );
            if( currentRetransmissionDelays_.size( ) != static_cast< unsigned int >( numberOfLinkEnds_ - 2 ) )
            {
                throw std::runtime_error( "Error when calculating n-way range, retransmission delay vector size is inconsistent" );
            }
        }
        else
        {
            for( int i = 0; i < numberOfLinkEnds_; i++ )
            {
                currentRetransmissionDelays_.push_back( 0.0 );
            }
        }

        // Retrieve index of link end where to start.
        int startLinkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndAssociatedWithTime, numberOfLinkEnds_ );
        int currentDownIndex = startLinkEndIndex;

        // Define 'current time'
        TimeType currentLinkEndStartTime = time;

        // Move 'backwards' from reference link end to transmitter.
        while( currentDownIndex > 0 )
        {
            currentLightTime = lightTimeCalculators_.at( currentDownIndex - 1 )->calculateLightTimeWithLinkEndsStates(
                        currentReceiverStateOutput, currentTransmitterStateOutput,
                        currentLinkEndStartTime, 1 );

            // Add link-end times/states for current leg.
            linkEndStates[ 2 * ( currentDownIndex - 1 ) + 1 ] = currentReceiverStateOutput.template cast< double >( );
            linkEndStates[ 2 * ( currentDownIndex - 1 ) ] = currentTransmitterStateOutput.template cast< double >( );
            linkEndTimes[ 2 * ( currentDownIndex - 1 ) + 1 ] = currentLinkEndStartTime;
            linkEndTimes[ 2 * ( currentDownIndex - 1 )] = currentLinkEndStartTime - currentLightTime;


            // If an additional leg is required, retrieve retransmission delay and update current time
            currentLinkEndStartTime -= currentLightTime;
            if( currentDownIndex > 1 )
            {
                currentLightTime += currentRetransmissionDelays_.at( currentDownIndex - 2 );
            }

            // Add computed light-time to total time and move to next leg
            totalLightTime += currentLightTime;
            currentDownIndex--;
        }

        int currentUpIndex = startLinkEndIndex;

        // If start is not at transmitter, compute and add retransmission delay.
        if( ( startLinkEndIndex != 0 ) && ( startLinkEndIndex != numberOfLinkEnds_ - 1 ) )
        {
            currentLinkEndStartTime = time + currentRetransmissionDelays_.at( startLinkEndIndex - 1 );
            totalLightTime += currentRetransmissionDelays_.at( startLinkEndIndex - 1 );
        }

        // Move 'forwards' from reference link end to receiver.
        while( currentUpIndex < static_cast< int >( lightTimeCalculators_.size( ) ) )
        {
            currentLightTime = lightTimeCalculators_.at( currentUpIndex )->calculateLightTimeWithLinkEndsStates(
                        currentReceiverStateOutput, currentTransmitterStateOutput,
                        currentLinkEndStartTime, 0 );

            // Add link-end times/states for current leg.
            linkEndStates[ 2 * currentUpIndex + 1 ] = currentReceiverStateOutput.template cast< double >( );
            linkEndStates[ 2 * currentUpIndex ] = currentTransmitterStateOutput.template cast< double >( );
            linkEndTimes[ 2 * currentUpIndex + 1 ] = currentLinkEndStartTime + currentLightTime;
            linkEndTimes[ 2 * currentUpIndex ] = currentLinkEndStartTime;

            // If an additional leg is required, retrieve retransmission delay and update current time
            currentLinkEndStartTime += currentLightTime;
            if( currentUpIndex < static_cast< int >( lightTimeCalculators_.size( ) ) - 1 )
            {
                currentLightTime += currentRetransmissionDelays_.at( currentUpIndex );
            }

            // Add computed light-time to total time and move to next leg
            totalLightTime += currentLightTime;
            currentUpIndex++;
        }

        // Return total range observation.
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >(
                     ) << totalLightTime * physical_constants::getSpeedOfLight< ObservationScalarType >( ) ).finished( );
    }

    std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > getLightTimeCalculators( )
    {
        return lightTimeCalculators_;
    }

private:

    //! List of objects to compute the light-times for each leg of the n-way range.
    /*!
     *  List of objects to compute the light-times (including any corrections w.r.t. Euclidean case)  for each leg of the
     *  n-way range.  First entry starts at transmitter; last entry is to receiver.
     */
    std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >
    lightTimeCalculators_;

    //!  Function that returns the list of retransmission delays as a function of observation time.
    /*!
     *  Function that returns the list of retransmission delays as a function of observation time.
     *  The retransmission delays represent the time difference between the reception of a singal by an intermediate link end,
     *  and the retransmission to the subsequent link end. By default, this function is empty, in which case no retransmission
     *  delays are used.
     */
    std::function< std::vector< double >( const double ) > retransmissionDelays_;

    //! List of retransmission delays, as computed by last call to computeIdealObservationsWithLinkEndData.
    std::vector< double > currentRetransmissionDelays_;

    //! Number of links in n-way observation
    int numberOfLinks_;

    //! Number of link ends in n-way observation (=number of links + 1)
    int numberOfLinkEnds_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_NWAYRANGEOBSERVATIONMODEL_H
