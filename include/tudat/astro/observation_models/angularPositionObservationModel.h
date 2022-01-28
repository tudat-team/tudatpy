/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ANGULARPOSITIONOBSERVATIONMODEL_H
#define TUDAT_ANGULARPOSITIONOBSERVATIONMODEL_H

#include <map>
#include <Eigen/Core>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/observationModel.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating angular position (right ascension/declination) observables.
/*!
 *  Class for simulating angular position (right ascension/declination), using light-time (with light-time corrections)
 *  to determine the states of the link ends (source and receiver).
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class AngularPositionObservationModel: public ObservationModel< 2, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > PositionType;

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between source and receiver
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    AngularPositionObservationModel(
            const LinkEnds linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
            const std::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = nullptr ):
        ObservationModel< 2, ObservationScalarType, TimeType >( angular_position, linkEnds, observationBiasCalculator ),
        lightTimeCalculator_( lightTimeCalculator ) { }

    //! Destructor
    ~AngularPositionObservationModel( ){ }

    //! Function to compute ideal angular position observation at given time.
    /*!
     *  This function compute ideal angular position observation at a given time. The time argument can be either the
     *  reception or transmission time (defined by linkEndAssociatedWithTime input).
     *  Note that this observable does include e.g. light-time corrections, which represent physically true corrections.
     *  It does not include e.g. system-dependent measurement.
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated angular position observable values.
     */
    Eigen::Matrix< ObservationScalarType, 2, 1 > computeIdealObservationsWithLinkEndData(
                    const TimeType time,
                    const LinkEndType linkEndAssociatedWithTime,
                    std::vector< double >& linkEndTimes,
                    std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )

    {
        // Check link end associated with input time and compute observable
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
            throw std::runtime_error( "Error when calculating angular position observation, link end is not transmitter or receiver" );
        }

        Eigen::Matrix< ObservationScalarType, 6, 1 > receiverState;
        Eigen::Matrix< ObservationScalarType, 6, 1 > transmitterState;

        // Compute light-time and receiver/transmitter states.
        ObservationScalarType lightTime = lightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, time, isTimeAtReception );

        // Compute spherical relative position
        Eigen::Matrix< ObservationScalarType, 3, 1 > sphericalRelativeCoordinates =
                coordinate_conversions::convertCartesianToSpherical< ObservationScalarType >(
                    transmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 ) ).
                template cast< ObservationScalarType >( );

        // Set link end times and states.
        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndStates.push_back( transmitterState.template cast< double >( ) );
        linkEndStates.push_back( receiverState.template cast< double >( ) );

        if( isTimeAtReception )
        {
            linkEndTimes.push_back( static_cast< double >( time - lightTime ) );
            linkEndTimes.push_back( static_cast< double >( time ) );
        }
        else
        {
            linkEndTimes.push_back( static_cast< double >( time ) );
            linkEndTimes.push_back( static_cast< double >( time + lightTime ) );
        }

        // Return observable
        return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << sphericalRelativeCoordinates.z( ),
                 mathematical_constants::PI / 2.0 - sphericalRelativeCoordinates.y( ) ).finished( );
    }

    //! Function to get the object to calculate light time.
    /*!
     * Function to get the object to calculate light time.
     * \return Object to calculate light time.
     */    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return lightTimeCalculator_;
    }

private:

    //! Object to calculate light time.
    /*!
     *  Object to calculate light time, including possible corrections from troposphere, relativistic corrections, etc.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator_;
};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_ANGULARPOSITIONOBSERVATIONMODEL_H
