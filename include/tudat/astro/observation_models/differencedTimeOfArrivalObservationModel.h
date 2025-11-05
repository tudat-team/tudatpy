/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIFFERENCEDTIMEOFARRIVALOBSERVATIONMODEL_H
#define TUDAT_DIFFERENCEDTIMEOFARRIVALOBSERVATIONMODEL_H

#include <map>

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{



template< typename ObservationScalarType = double, typename TimeType = double >
class OneWayDifferencedTimeOfArrivalObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;
    
    OneWayDifferencedTimeOfArrivalObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > firstReceiverLightTimeCalculator,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > secondReceiverLightTimeCalculator,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( differenced_time_of_arrival, linkEnds, observationBiasCalculator ),
        firstReceiverLightTimeCalculator_( firstReceiverLightTimeCalculator ), secondReceiverLightTimeCalculator_( secondReceiverLightTimeCalculator )
    {}

    //! Destructor
    ~OneWayDifferencedTimeOfArrivalObservationModel( ) {}
    
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetingsInput = nullptr )
    {
        ObservationScalarType lightTimeForFirstReceiver;
        ObservationScalarType lightTimeForSecondReceiver;

        linkEndTimes.resize( 3 );
        linkEndStates.resize( 3 );

        StateType transmitterStateForFirstLink, receiverStateForFirstLink, transmitterStateForSecondLink, receiverStateForSecondLink;
        if( linkEndAssociatedWithTime == receiver )
        {
            // Calculate reception time at ground station at the start and end of the count interval at reception time.
            linkEndTimes[ 1 ] = static_cast< double >( time );

            std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings;
            this->setFrequencyProperties( time, receiver, firstReceiverLightTimeCalculator_, ancilliarySetingsInput, ancilliarySetings );
            lightTimeForFirstReceiver = firstReceiverLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                    receiverStateForFirstLink, transmitterStateForFirstLink, linkEndTimes[ 1 ], 1, ancilliarySetings );
            
            linkEndTimes[ 0 ] = linkEndTimes[ 1 ] - static_cast< double >( lightTimeForFirstReceiver );

            this->setFrequencyProperties( time, receiver, secondReceiverLightTimeCalculator_, ancilliarySetingsInput, ancilliarySetings );
            lightTimeForSecondReceiver = secondReceiverLightTimeCalculator_->calculateLightTimeWithLinkEndsStates( 
                    receiverStateForSecondLink, transmitterStateForSecondLink, linkEndTimes[ 0 ], 0, ancilliarySetings );
            linkEndTimes[ 2 ] = linkEndTimes[ 0 ] + lightTimeForSecondReceiver;

        }
        else
        {
            throw std::runtime_error( "Error in differenced differenced time of arrival, reference link end not recognized" );
        }

        linkEndStates[ 0 ] = transmitterStateForFirstLink.template cast< double >( );
        linkEndStates[ 1 ] = receiverStateForFirstLink.template cast< double >( );
        linkEndStates[ 2 ] = receiverStateForSecondLink.template cast< double >( );

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << ( linkEndTimes[ 1 ] - linkEndTimes[ 2 ] ) ).finished( );
    }

    //! Light time calculator to compute light time at the beginning of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getFirstReceiverLightTimeCalculator( )
    {
        return firstReceiverLightTimeCalculator_;
    }

    //! Light time calculator to compute light time at the end of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getSecondReceiverLightTimeCalculator( )
    {
        return secondReceiverLightTimeCalculator_;
    }

private:
    //! Light time calculator to compute light time at the beginning of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > firstReceiverLightTimeCalculator_;

    //! Light time calculator to compute light time at the end of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > secondReceiverLightTimeCalculator_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_DIFFERENCEDTIMEOFARRIVALOBSERVATIONMODEL_H
