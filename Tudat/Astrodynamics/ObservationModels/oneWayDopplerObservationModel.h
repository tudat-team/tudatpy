#ifndef ONEWAYDOPPLEROBSERVATIONMODEL_H
#define ONEWAYDOPPLEROBSERVATIONMODEL_H


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

template< typename ObservationScalarType >
ObservationScalarType calculateLineOfSightVelocityAsCFraction(
        const Eigen::Matrix< ObservationScalarType, 3, 1 >& lineOfSightUnitVector,
        const Eigen::Matrix< ObservationScalarType, 3, 1 >& velocityVector )
{
    return lineOfSightUnitVector.dot( velocityVector ) / physical_constants::getSpeedOfLight< ObservationScalarType >( );
}

template< typename ObservationScalarType = double, typename TimeType = double >
class OneWayDopplerObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
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
    OneWayDopplerObservationModel(
            const boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL ):
        ObservationModel< 1, ObservationScalarType, TimeType >( oneWayDoppler, observationBiasCalculator ),
        lightTimeCalculator_( lightTimeCalculator )
    {
        one_ = mathematical_constants::getFloatingInteger< ObservationScalarType >( 1 );
    }

    ~OneWayDopplerObservationModel( ){ }

    //! Function to compute ideal one-way Doppler observation at given time.
    /*!
     *  This function compute ideal the one-way observation at a given time. The time argument can be either the reception
     *  or transmission time (defined by linkEndAssociatedWithTime input) Note that this observable does include e.g.
     *  light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \return Calculated observed one-way Doppler value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )

    {
        ObservationScalarType lightTime;
        StateType receiverState, transmitterState;

        bool isTimeAtReception = -1;

        switch( linkEndAssociatedWithTime )
        {
        case receiver:
            isTimeAtReception = true;
            break;
        case transmitter:
            isTimeAtReception = false;
            break;
        default:
            throw std::runtime_error(
                        "Error when calculating one way Doppler observation, link end is not transmitter or receiver" );
        }

        lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, time, true );

        PositionType relativePostion = ( ( receiverState - transmitterState ).segment( 0, 3 ) ).normalized( );

        return ( Eigen::Matrix<  ObservationScalarType, 1, 1  >( )
                <<( one_ - calculateLineOfSightVelocityAsCFraction< ObservationScalarType >( relativePostion, receiverState.segment( 3, 3 ) ) ) /
                ( one_ - calculateLineOfSightVelocityAsCFraction< ObservationScalarType >( relativePostion, transmitterState.segment( 3, 3 ) ) ) ).finished( );
    }

    //! Function to compute one-way Doppler observable without any corrections.
    /*!
     *  Function to compute one-way Doppler  observable without any corrections, i.e. the true physical Doppler as computed
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
     *  \return Ideal one-way Doppler observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        ObservationScalarType lightTime;
        StateType receiverState, transmitterState;
        TimeType transmissionTime, receptionTime;

        bool fixTransmissionTime;

        switch( linkEndAssociatedWithTime )
        {
        case receiver:
            lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState, transmitterState, time, true );
            transmissionTime = time - lightTime;
            receptionTime = time;
            fixTransmissionTime = 0;
            break;

        case transmitter:
            lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState, transmitterState, time, false );
            transmissionTime = time;
            receptionTime = time + lightTime;
            fixTransmissionTime = 1;
            break;
        default:
            throw std::runtime_error(
                        "Error when calculating one way Doppler observation, link end is not transmitter or receiver" );
        }

        linkEndTimes.push_back( transmissionTime );
        linkEndTimes.push_back( receptionTime );

        linkEndStates.push_back( transmitterState );
        linkEndStates.push_back( receiverState );

        PositionType relativePostion = ( ( receiverState - transmitterState ).segment( 0, 3 ) ).normalized( );
        relativePostion /= relativePostion.norm( );

        return ( Eigen::Matrix<  ObservationScalarType, 1, 1  >( )
                <<( one_ - calculateLineOfSightVelocityAsCFraction< ObservationScalarType >( relativePostion, receiverState.segment( 3, 3 ) ) ) /
                ( one_ - calculateLineOfSightVelocityAsCFraction< ObservationScalarType >( relativePostion, transmitterState.segment( 3, 3 ) ) ) ).finished( );

    }

    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }


private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator_;

    ObservationScalarType one_;

};

}

}

#endif // ONEWAYDOPPLEROBSERVATIONMODEL_H
