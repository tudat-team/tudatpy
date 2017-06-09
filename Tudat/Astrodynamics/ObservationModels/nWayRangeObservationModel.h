/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h>

namespace tudat
{

namespace observation_models
{

//! Class for simulating one-way range observables.
/*!
 *  Class for simulating one-way range, based on light-time and light-time corrections.
 *  The one-way range is defined as the light time multiplied by speed of light.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
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
     *  \param lightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    NWayRangeObservationModel(
            const std::vector< boost::shared_ptr< observation_models::LightTimeCalculator
            < ObservationScalarType, TimeType > > > lightTimeCalculators,
            const boost::function< std::vector< double >( ) > retransmissionDelays =
            boost::function< std::vector< double >( ) >( ),
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_range, observationBiasCalculator ),
        lightTimeCalculators_( lightTimeCalculators ), retransmissionDelays_( retransmissionDelays )
    {
        numberOfLinks_ = lightTimeCalculators_.size( );
        numberOfLinkEnds_ = numberOfLinks_ + 1;
    }

    //! Destructor
    ~NWayRangeObservationModel( ){ }

    //! Function to compute ideal one-way range observation at given time.
    /*!
     *  This function compute ideal the one-way observation at a given time. The time argument can be either the reception
     *  or transmission time (defined by linkEndAssociatedWithTime input) Note that this observable does include e.g.
     *  light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \return Calculated observed one-way range value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )

    {
        return computeIdealObservationsWithLinkEndData( time, linkEndAssociatedWithTime, linkEndTimes_, linkEndStates_ );
    }

    //! Function to compute one-way range observable without any corrections.
    /*!
     *  Function to compute one-way range  observable without any corrections, i.e. the true physical range as computed
     *  from the defined link ends. Note that this observable does include light-time
     *  corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way range observable.
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
        if( !retransmissionDelays_.empty( ) )
        {
            currentRetransmissionDelays_ = retransmissionDelays_( );
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


            linkEndStates[ 2 * currentUpIndex + 1 ] = currentReceiverStateOutput.template cast< double >( );
            linkEndStates[ 2 * currentUpIndex ] = currentTransmitterStateOutput.template cast< double >( );

            linkEndTimes[ 2 * currentUpIndex + 1 ] = currentLinkEndStartTime + currentLightTime;
            linkEndTimes[ 2 * currentUpIndex ] = currentLinkEndStartTime;

            if( currentUpIndex < static_cast< int >( lightTimeCalculators_.size( ) ) - 1 )
            {
                currentLightTime += currentRetransmissionDelays_.at( currentUpIndex );
            }

            currentLinkEndStartTime += currentLightTime;

            totalLightTime += currentLightTime;
            currentUpIndex++;
        }

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >(
                     ) <<totalLightTime * physical_constants::getSpeedOfLight< ObservationScalarType >( ) ).finished( );
    }

    std::vector< boost::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > getLightTimeCalculators( )
    {
        return lightTimeCalculators_;
    }

private:

    std::vector< boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >
    lightTimeCalculators_;

    boost::function< std::vector< double >( ) > retransmissionDelays_;

    std::vector< double > currentRetransmissionDelays_;

    int numberOfLinks_;

    int numberOfLinkEnds_;

    //! Pre-declared vector of link end times, used for computeIdealObservations function
    std::vector< double > linkEndTimes_;

    //! Pre-declared vector of link end states, used for computeIdealObservations function
    std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_NWAYRANGEOBSERVATIONMODEL_H
