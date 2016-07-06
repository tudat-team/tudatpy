/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef ANGULARPOSITIONOBSERVATIONMODEL_H
#define ANGULARPOSITIONOBSERVATIONMODEL_H

#include <map>
#include <string>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/ObservationModels/createLightTimeCalculator.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class AngularPositionObservationModel: public ObservationModel< 2, ObservationScalarType, TimeType, StateScalarType >
{
public:

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< StateScalarType, 6, 1 > PositionType;


    AngularPositionObservationModel(
            const boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > lightTimeCalculator,
            const boost::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = NULL ):
        ObservationModel< 2, ObservationScalarType, TimeType, StateScalarType >( angular_position, observationBiasCalculator ),
        lightTimeCalculator_( lightTimeCalculator ) { }


    Eigen::Matrix< ObservationScalarType, 2, 1 > computeUnbiasedObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< TimeType >& linkEndTimes,
                    std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const

    {
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
            isTimeAtReception = -1;
            std::cerr<<"Error when calculating angular position observation, link end is not transmitter or receiver"<<std::endl;
        }

        Eigen::Matrix< StateScalarType, 6, 1 > receiverState;
        Eigen::Matrix< StateScalarType, 6, 1 > transmitterState;

        ObservationScalarType lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, time, isTimeAtReception );

        Eigen::Matrix< StateScalarType, 3, 1 > vectorToTransmitter = transmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 );
        Eigen::Matrix< ObservationScalarType, 3, 1 > sphericalRelativeCoordinates =
                coordinate_conversions::convertCartesianToSpherical< StateScalarType >( vectorToTransmitter ).
                template cast< ObservationScalarType >( );

        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndStates.push_back( transmitterState );
        linkEndStates.push_back( receiverState );

        if( isTimeAtReception )
        {
            linkEndTimes.push_back( time - lightTime );
            linkEndTimes.push_back( time );
        }
        else
        {
            linkEndTimes.push_back( time );
            linkEndTimes.push_back( time + lightTime );
        }

        return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << sphericalRelativeCoordinates.z( ),
                 mathematical_constants::PI / 2.0 - sphericalRelativeCoordinates.y( ) ).finished( );
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

#endif // ANGULARPOSITIONOBSERVATIONMODEL_H
