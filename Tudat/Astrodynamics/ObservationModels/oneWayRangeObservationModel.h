#ifndef SIMPLERANGEOBSERVATIONMODEL_H
#define SIMPLERANGEOBSERVATIONMODEL_H

#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/createLightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating one-way range observables.
/*!
 *  Class for simulating one-way range, based on light-time and light-time corrections.
 *  The one-way range is defined as the observed 1-way range, i.e. light time multiplied by speed of
 *  light. The user may add observation biases to model system-dependent deviations between measured
 *  and true osbervation.
 */
template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class OneWayRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType >
{
public:

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< StateScalarType, 6, 1 > PositionType;


    //! Constructor.
    /*!
     *  Constructor, takes data defining the states of the linke ends and required corrections.
     *  \param transmitterCompleteEphemeris State function for the transmitter.
     *  \param receiverCompleteEphemeris State function for the receiver.
     *  \param lightTimeCorrections List of settings for light-time corrections (default is none).
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     * observable, i.e. deviations from the physically true observable (default none).
     */
    OneWayRangeObservationModel(
            const boost::function< StateType( const TimeType ) > transmitterCompleteEphemeris,
            const boost::function< StateType( const TimeType ) > receiverCompleteEphemeris,
            const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections =
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
            const boost::shared_ptr< ObservationBiasInterface > observationBiasCalculator = NULL ):
        ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType >( oneWayRange, observationBiasCalculator )
    {
        lightTimeCalculator_ = createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
                    transmitterCompleteEphemeris, receiverCompleteEphemeris, bodyMap, lightTimeCorrections, transmitter, receiver );
    }

    //! Destructor
    ~OneWayRangeObservationModel( ){ }

    //! Function to compute one-way range observation at given time.
    /*!
     *  This function computes the one-way observation at a given time.
     *  The time argument can be either the reception or transmission time.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which current time is measured, i.e. reference
     *  link end for observable.
     *  \return Calculated observed one-way range value.
     */
    ObservationScalarType computeObservation( const TimeType time,
                                              const LinkEndType linkEndAssociatedWithTime ) const

    {
        // Check link end associated with input time.
        bool isTimeAtReception;
        if( linkEndAssociatedWithTime == receiver )
        {
            isTimeAtReception = 1;
        }
        else if( linkEndAssociatedWithTime == transmitter )
        {
            isTimeAtReception = 0;
        }
        else
        {
            std::cerr<<"Error when calculating one way range observation, link end is not transmitter or receiver"<<std::endl;
        }

        // Calculate light-time and multiply by speed of light in vacuum.
        return computeLightTime( time, isTimeAtReception ) *
                physical_constants::getSpeedOfLight< ObservationScalarType >( );
    }

    ObservationScalarType computeObservationAndFullPrecisionLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< TimeType >& linkEndTimes,
            std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const
    {
        ObservationScalarType observation;
        StateType receiverState, transmitterState;
        TimeType transmissionTime, receptionTime;
        switch( linkEndAssociatedWithTime )
        {
        case receiver:
            observation = computeLightTimeFromReceptionWithLinkEndsStates(
                        time, receiverState, transmitterState );
            transmissionTime = time - observation;
            receptionTime = time;
            break;

        case transmitter:
            observation = computeLightTimeFromTransmissionWithLinkEndsStates(
                        time, receiverState, transmitterState );
            transmissionTime = time;
            receptionTime = time + observation;
            break;
        }



        observation *= physical_constants::getSpeedOfLight< ObservationScalarType >( );

        linkEndTimes.clear( );
        linkEndStates.clear( );

        linkEndTimes.push_back( transmissionTime );
        linkEndTimes.push_back( receptionTime );

        linkEndStates.push_back( transmitterState );
        linkEndStates.push_back( receiverState );

        return observation;
    }

    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > lightTimeCalculator_;

};

}

}

#endif // SIMPLERANGEOBSERVATIONMODEL_H
