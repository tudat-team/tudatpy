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
     *  \param observationBiasCalculator Object for calculating (system-dependent) errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    NWayRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::vector< std::shared_ptr< observation_models::LightTimeCalculator
                < ObservationScalarType, TimeType > > > lightTimeCalculators,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
                = std::make_shared< MultiLegLightTimeConvergenceCriteria >( ) ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_range, linkEnds, observationBiasCalculator )
    {
        multiLegLightTimeCalculator_ = std::make_shared< observation_models::MultiLegLightTimeCalculator<
                ObservationScalarType, TimeType > >( lightTimeCalculators, lightTimeConvergenceCriteria );
    }

    NWayRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::MultiLegLightTimeCalculator
                < ObservationScalarType, TimeType > > multiLegLightTimeCalculator,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( n_way_range, linkEnds, observationBiasCalculator ),
        multiLegLightTimeCalculator_( multiLegLightTimeCalculator )
    { }

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
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings = nullptr  )
    {

        ObservationScalarType totalLightTime = multiLegLightTimeCalculator_->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancilliarySetings );

        // Return total range observation.
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >(
                     ) << totalLightTime * physical_constants::getSpeedOfLight< ObservationScalarType >( ) ).finished( );
    }

    std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > getLightTimeCalculators( )
    {
        return multiLegLightTimeCalculator_->getLightTimeCalculators( );
    }

    std::shared_ptr< MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > getMultiLegLightTimeCalculator( )
    {
        return multiLegLightTimeCalculator_;
    }

private:

    // Object that iteratively computes the light time of multiple legs
    std::shared_ptr< MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > multiLegLightTimeCalculator_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_NWAYRANGEOBSERVATIONMODEL_H
