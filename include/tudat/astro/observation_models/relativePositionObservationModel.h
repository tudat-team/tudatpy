/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVEPOSITIONOBSERVATIONMODEL_H
#define TUDAT_RELATIVEPOSITIONOBSERVATIONMODEL_H


#include <functional>

#include "tudat/astro/ephemerides/ephemeris.h"

#include "tudat/astro/observation_models/observationModel.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating observations of three-dimensional position.
/*!
 *  Class for simulating observations of three-dimensional position. This observable is typically not realized in
 *  practice, but its use can be very valuable in simulation studies
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class RelativePositionObservationModel: public ObservationModel< 3, ObservationScalarType, TimeType >
{
public:

    //! Constructor.
    /*!
     *  Constructor,
     *  \param stateFunction Function that returns the Cartesian state of the observed body as a function of time.
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable (default none).
     */
    RelativePositionObservationModel(
            const LinkEnds& linkEnds,
            const std::function<  Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType& ) > stateFunctionTarget,
            const std::function<  Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType& ) > stateFunctionObserver,
            const std::shared_ptr< ObservationBias< 3 > > observationBiasCalculator = nullptr ):
        ObservationModel< 3, ObservationScalarType, TimeType >(
            relative_position_observable, linkEnds, observationBiasCalculator ), stateFunctionTarget_( stateFunctionTarget ), stateFunctionObserver_( stateFunctionObserver ){ }

    //! Destructor
    ~RelativePositionObservationModel( ) { }

    //! Function to compute ideal position observation at given time.
    /*!
     *  This function computes the ideal position observation at a given time (without biases).
     *  The times and states of the link ends are also returned in full precision (determined by class template
     *  arguments). These states and times are returned by reference.
     * \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid (must be observed_body for this derived class)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal position observable.
     */
    Eigen::Matrix< ObservationScalarType, 3, 1 > computeIdealObservationsWithLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< double >& linkEndTimes,
                std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
                const std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySetings = nullptr )
    {
        // Check link end
        if( ( linkEndAssociatedWithTime != observed_body ) && ( linkEndAssociatedWithTime != observer ) )
        {
            throw std::runtime_error(
                        "Error when computing position observable, associated link end must be observed_body or observer." );
        }

        currentStateTarget_ = stateFunctionTarget_( time );
        currentStateObserver_ = stateFunctionObserver_( time );

        // Set link end times and states.
        linkEndTimes.clear( );
        linkEndTimes.push_back( static_cast< double >( time ) );
        linkEndTimes.push_back( static_cast< double >( time ) );

        linkEndStates.clear( );
        linkEndStates.push_back( currentStateTarget_.template cast< double >( ) );
        linkEndStates.push_back( currentStateObserver_.template cast< double >( ) );

        // Retrieve position
        return ( currentStateTarget_.segment( 0, 3 ) - currentStateObserver_.segment( 0, 3 ) );
    }


private:

    //! Function that returns the Cartesian state of the observed body as a function of time.
    std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType& ) > stateFunctionTarget_;

    //! Function that returns the Cartesian state of the observer body as a function of time.
    std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType& ) > stateFunctionObserver_;

    Eigen::Matrix< ObservationScalarType, 6, 1 > currentStateTarget_;

    Eigen::Matrix< ObservationScalarType, 6, 1 > currentStateObserver_;
};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_RELATIVEPOSITIONOBSERVATIONMODEL_H
