/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONEWAYRANGEOBSERVATIONMODEL_H
#define TUDAT_ONEWAYRANGEOBSERVATIONMODEL_H

#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/ObservationModels/createLightTimeCalculator.h"
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
            const boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > lightTimeCalculator,
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL ):
        ObservationModel< 1, ObservationScalarType, TimeType, StateScalarType >( oneWayRange, observationBiasCalculator ),
      lightTimeCalculator_( lightTimeCalculator ){ }

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
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeUnbiasedObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime ) const

    {
        // Check link end associated with input time.
        bool isTimeAtReception = -1;
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
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) <<
                 lightTimeCalculator_->calculateLightTime( time, isTimeAtReception ) *
                 physical_constants::getSpeedOfLight< ObservationScalarType >( ) ).finished( );
    }

    Eigen::Matrix< ObservationScalarType, 1, 1 > computeUnbiasedObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< TimeType >& linkEndTimes,
                    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const
    {
        ObservationScalarType observation = TUDAT_NAN;
        StateType receiverState, transmitterState;
        TimeType transmissionTime, receptionTime;
        switch( linkEndAssociatedWithTime )
        {
        case receiver:
            observation = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState, transmitterState, time, 1 );
            transmissionTime = time - observation;
            receptionTime = time;
            break;

        case transmitter:
            observation = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                        receiverState, transmitterState, time, 0 );
            transmissionTime = time;
            receptionTime = time + observation;
            break;
        default:
            std::string errorMessage = "Error, cannot have link end type: " +
                    boost::lexical_cast< std::string >( linkEndAssociatedWithTime ) + "for one-way range";
            throw std::runtime_error( errorMessage );
        }

        observation *= physical_constants::getSpeedOfLight< ObservationScalarType >( );


        linkEndTimes.push_back( transmissionTime );
        linkEndTimes.push_back( receptionTime );

        linkEndStates.push_back( transmitterState);
        linkEndStates.push_back( receiverState );

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << observation ).finished( );
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

#endif // TUDAT_ONEWAYRANGEOBSERVATIONMODEL_H
