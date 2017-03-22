#ifndef TUDAT_ONEWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H
#define TUDAT_ONEWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H


#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

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
     *  \param dopplerIntervalFunction Function returning the integration time of the observable as a function of the
     *  current observation time.
     */
    OneWayDifferencedRangeObservationModel(
            const boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
            arcStartLightTimeCalculator,
            const boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
            arcEndLightTimeCalculator,
            boost::function< double( const double ) > dopplerIntervalFunction,
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL ):
        ObservationModel< 1, ObservationScalarType, TimeType >( one_way_differenced_range, observationBiasCalculator ),
        arcStartLightTimeCalculator_( arcStartLightTimeCalculator ), arcEndLightTimeCalculator_( arcEndLightTimeCalculator ),
        dopplerIntervalFunction_( dopplerIntervalFunction )
    {

    }

    //! Destructor
    ~OneWayDifferencedRangeObservationModel( ){ }

    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )

    {
        ObservationScalarType lightTimeAtStartInterval;
        ObservationScalarType lightTimeAtEndInterval;
        TimeType currentDopplerCountInterval = dopplerIntervalFunction_( time );

        if ( linkEndAssociatedWithTime == receiver )
        {
            // Calculate light times at the start of the reception interval
            lightTimeAtStartInterval = arcStartLightTimeCalculator_->calculateLightTime(
                        static_cast< double >( time ) - currentDopplerCountInterval, 1 );

            // Calculate light times at the end of the reception interval
            lightTimeAtEndInterval = arcEndLightTimeCalculator_->calculateLightTime(
                        static_cast< double >( time ) , 1 );

        }
        else if ( linkEndAssociatedWithTime == transmitter )
        {

            // Calculate light times at the start of the reception interval
            lightTimeAtEndInterval = arcEndLightTimeCalculator_->calculateLightTime(
                        static_cast< double >( time ) - currentDopplerCountInterval, 0 );

            // Calculate light times at the end of the reception interval
            lightTimeAtStartInterval = arcStartLightTimeCalculator_->calculateLightTime(
                        static_cast< double >( time ), 0 );

        }
        else
        {
            throw std::runtime_error( "Error in differenced range rate observation model, reference link end not recognized" );
        }

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << ( lightTimeAtEndInterval - lightTimeAtStartInterval ) *
                 physical_constants::getSpeedOfLight< ObservationScalarType >( ) /
                 static_cast< ObservationScalarType >( currentDopplerCountInterval ) ).finished( );
    }

    //! Function to compute one-way Doppler observable without any corrections.
    /*!
     *  Function to compute one-way Doppler  observable without any corrections, i.e. the true physical Doppler as computed
     *  from the defined link ends. It does not include system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in double precision. These states and times are returned by
     *  reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal one-way Doppler observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        ObservationScalarType lightTimeAtStartInterval;
        ObservationScalarType lightTimeAtEndInterval;
        TimeType currentDopplerCountInterval = dopplerIntervalFunction_( time );

        linkEndTimes.resize( 4 );
        linkEndStates.resize( 4 );

        StateType transmitterStateAtArcStart, receiverStateAtArcStart, transmitterStateAtArcEnd, receiverStateAtArcEnd;
        if ( linkEndAssociatedWithTime == receiver )
        {
            //Calculate reception time at ground station at the start and end of the count interval at reception time.
            linkEndTimes[ 1 ] = static_cast< double >( time ) - currentDopplerCountInterval;
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
            linkEndTimes[ 0 ] = static_cast< double >( time ) - currentDopplerCountInterval;
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
                 static_cast< ObservationScalarType >( currentDopplerCountInterval ) ).finished( );
    }


private:

    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
    arcStartLightTimeCalculator_;

    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
    arcEndLightTimeCalculator_;

    boost::function< double( const double ) > dopplerIntervalFunction_;

};

}

}

#endif // TUDAT_ONEWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H
