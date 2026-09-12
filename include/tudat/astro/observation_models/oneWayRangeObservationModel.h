/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <functional>

#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/transmissionFrequencyInterface.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating one-way range observables.
/*!
 *  Class for simulating one-way range, based on light-time and light-time corrections.
 *  The one-way range is defined as the light time multiplied by speed of light.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class OneWayRangeObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
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
    OneWayRangeObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const basic_astrodynamics::TimeScales scaleForTimeDifference = basic_astrodynamics::tdb_scale,
            const std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > groundStationStates =
                    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >( ) ):
        ObservationModel< 1, ObservationScalarType, TimeType >(
                one_way_range,
                linkEnds,
                observationBiasCalculator,
                std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >{ std::make_shared<
                        FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                        std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{
                                lightTimeCalculator },
                        std::make_shared< LightTimeConvergenceCriteria >( ),
                        false ) } ),
        scaleForTimeDifference_( scaleForTimeDifference ), stationStates_( groundStationStates )
    {}

    //! Destructor
    ~OneWayRangeObservationModel( ) {}

    //! Function to compute one-way range observable without any corrections.
    /*!
     *  Function to compute one-way range  observable without any corrections, i.e. the true physical range as computed
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
     *  \return Ideal one-way range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetings = nullptr ) override
    {
        std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetingsToUse;
        std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator =
                getLightTimeCalculator( );
        this->setFrequencyProperties( time, linkEndAssociatedWithTime, lightTimeCalculator, ancillarySetings, ancillarySetingsToUse );

        ObservationScalarType observation = this->getFullLinkLightTimeCalculatorFromBase( )->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySetingsToUse );

        if( scaleForTimeDifference_ != basic_astrodynamics::tdb_scale )
        {
            Eigen::Vector3d nominalReceivingStationState = ( stationStates_.count( receiver ) == 0 )
                    ? Eigen::Vector3d::Zero( )
                    : stationStates_.at( receiver )->getNominalCartesianPosition( );

            Eigen::Vector3d nominalTransmittingStationState = ( stationStates_.count( transmitter ) == 0 )
                    ? Eigen::Vector3d::Zero( )
                    : stationStates_.at( transmitter )->getNominalCartesianPosition( );

            TimeType transmissionTimeDifference = this->timeScaleConverter_->getCurrentTimeDifference(
                    basic_astrodynamics::tdb_scale, scaleForTimeDifference_, linkEndTimes.front( ), nominalTransmittingStationState );
            TimeType receptionTimeDifference = this->timeScaleConverter_->getCurrentTimeDifference(
                    basic_astrodynamics::tdb_scale, scaleForTimeDifference_, linkEndTimes.back( ), nominalReceivingStationState );
            observation += convertIndependentVariableToScalar< ObservationScalarType >( receptionTimeDifference );
            observation -= convertIndependentVariableToScalar< ObservationScalarType >( transmissionTimeDifference );
        }

        // Convert light time to range.
        observation *= physical_constants::getSpeedOfLight< ObservationScalarType >( );

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << observation ).finished( );
    }

    //! Function to get the object to calculate light time.
    /*!
     * Function to get the object to calculate light time.
     * \return Object to calculate light time.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculator( )
    {
        return this->getSingleLegLightTimeCalculator( );
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > >
    getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( ) } } };
    }

private:
    basic_astrodynamics::TimeScales scaleForTimeDifference_;

    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_ONEWAYRANGEOBSERVATIONMODEL_H
