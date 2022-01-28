/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONEWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H
#define TUDAT_ONEWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H

#include <map>

#include <functional>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating one-way differenced range (e.g. closed-loop Doppler) observable
/*!
 *  Class for simulating one-way differenced range (e.g. closed-loop Doppler) observable. The observable is obtained by
 *  subtracting the range at two time intervals, and dividing by the time difference. It represents the time-averages value
 *  of the range-rate over this integration time.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class OneWayDifferencedRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;


    //! Constructor.
    /*!
     *  Constructor,
     *  \param arcStartLightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param arcEndLightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     *  \param integrationTimeFunction Function returning the integration time of the observable as a function of the
     *  current observation time.
     */
    OneWayDifferencedRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
            arcStartLightTimeCalculator,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
            arcEndLightTimeCalculator,
            std::function< double( const double ) > integrationTimeFunction,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( one_way_differenced_range, linkEnds, observationBiasCalculator ),
        arcStartLightTimeCalculator_( arcStartLightTimeCalculator ), arcEndLightTimeCalculator_( arcEndLightTimeCalculator ),
        integrationTimeFunction_( integrationTimeFunction )
    {

    }

    //! Destructor
    ~OneWayDifferencedRangeObservationModel( ){ }

    //! Function to compute one-way differenced range observable without any corrections.
    /*!
     *  Function to compute one-way differenced range  observable without any corrections, i.e. the true physical differenced
     *  range as computed from the defined link ends. It does not include system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in double precision. These states and times are returned by
     *  reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way differenced range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        ObservationScalarType lightTimeAtStartInterval;
        ObservationScalarType lightTimeAtEndInterval;
        TimeType currentIntegrationTime = integrationTimeFunction_( time );

        linkEndTimes.resize( 4 );
        linkEndStates.resize( 4 );

        StateType transmitterStateAtArcStart, receiverStateAtArcStart, transmitterStateAtArcEnd, receiverStateAtArcEnd;
        if ( linkEndAssociatedWithTime == receiver )
        {
            //Calculate reception time at ground station at the start and end of the count interval at reception time.
            linkEndTimes[ 1 ] = static_cast< double >( time ) - currentIntegrationTime;
            linkEndTimes[ 3 ] = static_cast< double >( time );

            // Calculate light times at the start of the reception interval
            lightTimeAtStartInterval = arcStartLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                         receiverStateAtArcStart, transmitterStateAtArcStart, linkEndTimes[ 1 ] , 1 );

            // Calculate light times at the end of the reception interval
            lightTimeAtEndInterval = arcEndLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverStateAtArcEnd, transmitterStateAtArcEnd, linkEndTimes[ 3 ] , 1 );

            linkEndTimes[ 0 ] = linkEndTimes[ 1 ] - static_cast< double >( lightTimeAtStartInterval );
            linkEndTimes[ 2 ] = linkEndTimes[ 3 ] - static_cast< double >( lightTimeAtEndInterval );

        }
        else if ( linkEndAssociatedWithTime == transmitter )
        {
            //Calculate reception time at ground station at the start and end of the count interval at reception time.
            linkEndTimes[ 0 ] = static_cast< double >( time ) - currentIntegrationTime;
            linkEndTimes[ 2 ] = static_cast< double >( time );

            // Calculate light times at the start of the reception interval
            lightTimeAtEndInterval = arcEndLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverStateAtArcEnd, transmitterStateAtArcEnd, linkEndTimes[ 2 ], 0 );

            linkEndTimes[ 3 ] = linkEndTimes[ 2 ] + static_cast< double >( lightTimeAtEndInterval );

            // Calculate light times at the end of the reception interval
            lightTimeAtStartInterval = arcStartLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverStateAtArcStart, transmitterStateAtArcStart, linkEndTimes[ 0 ], 0 );

            linkEndTimes[ 1 ] = linkEndTimes[ 0 ] + static_cast< double >( lightTimeAtStartInterval );
        }
        else
        {
            throw std::runtime_error( "Error in differenced range rate observation model, reference link end not recognized" );
        }

        linkEndStates[ 0 ] = transmitterStateAtArcStart.template cast< double >( );
        linkEndStates[ 1 ] = receiverStateAtArcStart.template cast< double >( );
        linkEndStates[ 2 ] = transmitterStateAtArcEnd.template cast< double >( );
        linkEndStates[ 3 ] = receiverStateAtArcEnd.template cast< double >( );

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << ( lightTimeAtEndInterval - lightTimeAtStartInterval ) *
                 physical_constants::getSpeedOfLight< ObservationScalarType >( ) /
                 static_cast< ObservationScalarType >( currentIntegrationTime ) ).finished( );
    }

    //! Light time calculator to compute light time at the beginning of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getArcStartLightTimeCalculator( )
    {
        return arcStartLightTimeCalculator_;
    }

    //! Light time calculator to compute light time at the end of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getArcEndLightTimeCalculator( )
    {
        return arcEndLightTimeCalculator_;
    }

private:

    //! Light time calculator to compute light time at the beginning of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
    arcStartLightTimeCalculator_;

    //! Light time calculator to compute light time at the end of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
    arcEndLightTimeCalculator_;

    //! Function returning the integration time of the observable as a function of the current observation time.
    std::function< double( const double ) > integrationTimeFunction_;

};

}

}

#endif // TUDAT_ONEWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H
