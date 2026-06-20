/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONANGLEOBSERVATIONMODEL_H
#define TUDAT_POSITIONANGLEOBSERVATIONMODEL_H

#include <map>
#include <Eigen/Core>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/observationModel.h"

namespace tudat
{

namespace observation_models
{

inline double getPositionAngleScalingFactor( const observation_models::LinkEndType referenceLinkEnd,
                                             const std::vector< Eigen::Vector6d >& linkEndStates,
                                             const std::vector< double >& linkEndTimes,
                                             const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings,
                                             const bool isFirstPartial )
{
    return 1.0;
}

//! Class for simulating position angle observables.
/*!
 *  Class for simulating position angle observables, using light-time (with light-time corrections)
 *  to determine the states of the link ends (two transmitters and receiver).
 *  The position angle is measured from north through east, i.e. from the direction of the first
 *  transmitter towards the second transmitter.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class PositionAngleObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    static std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >
    createFullLinkLightTimeCalculators( const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                                lightTimeCalculatorFirstTransmitter,
                                        const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                                lightTimeCalculatorSecondTransmitter )
    {
        return std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >{
            std::make_shared< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                    std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{
                            lightTimeCalculatorFirstTransmitter },
                    std::make_shared< LightTimeConvergenceCriteria >( ),
                    false ),
            std::make_shared< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                    std::vector< std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{
                            lightTimeCalculatorSecondTransmitter },
                    std::make_shared< LightTimeConvergenceCriteria >( ),
                    false )
        };
    }

    //! Constructor.
    /*!
     *  Constructor,
     *  \param lightTimeCalculatorFirstTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between first transmitter and receiver
     *  \param lightTimeCalculatorSecondTransmitter Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  between second transmitter and receiver
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    PositionAngleObservationModel( const LinkEnds linkEnds,
                                   const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                           lightTimeCalculatorFirstTransmitter,
                                   const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                           lightTimeCalculatorSecondTransmitter,
                                   const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >(
                position_angle,
                linkEnds,
                observationBiasCalculator,
                createFullLinkLightTimeCalculators( lightTimeCalculatorFirstTransmitter, lightTimeCalculatorSecondTransmitter ) )
    {}

    //! Destructor
    ~PositionAngleObservationModel( ) {}

    //! Function to compute ideal position angle observation at given time, between two transmitters.
    /*!
     *  This function compute ideal position angle observation at a given time, between two transmitters. The time argument can
     * be either the reception or transmission time (defined by linkEndAssociatedWithTime input). Note that this observable does include
     * e.g. light-time corrections, which represent physically true corrections. It does not include e.g. system-dependent measurement. The
     * times and states of the link ends are also returned in full precision (determined by class template arguments). These states and
     * times are returned by reference.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated position angle observable value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettingsInput = nullptr ) override
    {
        // Check link end associated with input time and compute observable.
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error when calculating position angle observation, link end associated with time is not receiver." );
        }

        // Compute light-times and receiver/transmitters states.
        std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings;
        std::vector< double > firstLinkEndTimes;
        std::vector< double > secondLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > firstLinkEndStates;
        std::vector< Eigen::Matrix< double, 6, 1 > > secondLinkEndStates;
        std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorFirstTransmitter =
                getLightTimeCalculatorFirstTransmitter( );
        std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculatorSecondTransmitter =
                getLightTimeCalculatorSecondTransmitter( );
        this->setFrequencyProperties(
                time, linkEndAssociatedWithTime, lightTimeCalculatorFirstTransmitter, ancillarySettingsInput, ancillarySettings );
        this->getFullLinkLightTimeCalculatorFromBase( 0 )->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, firstLinkEndTimes, firstLinkEndStates, ancillarySettings );

        this->setFrequencyProperties(
                time, linkEndAssociatedWithTime, lightTimeCalculatorSecondTransmitter, ancillarySettingsInput, ancillarySettings );
        this->getFullLinkLightTimeCalculatorFromBase( 1 )->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, secondLinkEndTimes, secondLinkEndStates, ancillarySettings );

        Eigen::Matrix< ObservationScalarType, 6, 1 > receiverState = firstLinkEndStates.at( 1 ).template cast< ObservationScalarType >( );
        Eigen::Matrix< ObservationScalarType, 6, 1 > firstTransmitterState =
                firstLinkEndStates.at( 0 ).template cast< ObservationScalarType >( );
        Eigen::Matrix< ObservationScalarType, 6, 1 > secondTransmitterState =
                secondLinkEndStates.at( 0 ).template cast< ObservationScalarType >( );

        // Compute relative position vectors
        Eigen::Matrix< ObservationScalarType, 3, 1 > relativeStateTransmitter1 =
                firstTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 );
        Eigen::Matrix< ObservationScalarType, 3, 1 > relativeStateTransmitter2 =
                secondTransmitterState.segment( 0, 3 ) - receiverState.segment( 0, 3 );

        // Compute right ascension and declination for first transmitter.
        double rightAscensionFirstTransmitter = 2.0 *
                std::atan( relativeStateTransmitter1[ 1 ] /
                           ( std::sqrt( relativeStateTransmitter1[ 0 ] * relativeStateTransmitter1[ 0 ] +
                                        relativeStateTransmitter1[ 1 ] * relativeStateTransmitter1[ 1 ] ) +
                             relativeStateTransmitter1[ 0 ] ) );
        double declinationFirstTransmitter =
                mathematical_constants::PI / 2.0 - std::acos( relativeStateTransmitter1[ 2 ] / relativeStateTransmitter1.norm( ) );

        // Compute right ascension and declination for second transmitter.
        double rightAscensionSecondTransmitter = 2.0 *
                std::atan( relativeStateTransmitter2[ 1 ] /
                           ( std::sqrt( relativeStateTransmitter2[ 0 ] * relativeStateTransmitter2[ 0 ] +
                                        relativeStateTransmitter2[ 1 ] * relativeStateTransmitter2[ 1 ] ) +
                             relativeStateTransmitter2[ 0 ] ) );
        double declinationSecondTransmitter =
                mathematical_constants::PI / 2.0 - std::acos( relativeStateTransmitter2[ 2 ] / relativeStateTransmitter2.norm( ) );

        double deltaRA = rightAscensionSecondTransmitter - rightAscensionFirstTransmitter;

        // Compute position angle: atan2(sin(dRA) * cos(d2), cos(d1) * sin(d2) - sin(d1) * cos(d2) * cos(dRA))
        double positionAngle = std::atan2(
                std::sin( deltaRA ) * std::cos( declinationSecondTransmitter ),
                std::cos( declinationFirstTransmitter ) * std::sin( declinationSecondTransmitter ) -
                        std::sin( declinationFirstTransmitter ) * std::cos( declinationSecondTransmitter ) * std::cos( deltaRA ) );

        // Set link end times and states.
        linkEndTimes.clear( );
        linkEndStates.clear( );

        linkEndStates.push_back( firstLinkEndStates.at( 0 ) );
        linkEndStates.push_back( secondLinkEndStates.at( 0 ) );
        linkEndStates.push_back( firstLinkEndStates.at( 1 ) );

        linkEndTimes.push_back( firstLinkEndTimes.at( 0 ) );
        linkEndTimes.push_back( secondLinkEndTimes.at( 0 ) );
        linkEndTimes.push_back( firstLinkEndTimes.at( 1 ) );

        // Return observable
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << positionAngle ).finished( );
    }

    //! Function to get the object to calculate light time between first transmitter and receiver.
    /*!
     * Function to get the object to calculate light time between first transmitter and receiver.
     * \return Object to calculate light time between first transmitter and receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return this->getSingleLegLightTimeCalculator( 0, 0 );
    }

    //! Function to get the object to calculate light time between second transmitter and receiver.
    /*!
     * Function to get the object to calculate light time between second transmitter and receiver.
     * \return Object to calculate light time between second transmitter and receiver.
     */
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorSecondTransmitter( )
    {
        return this->getSingleLegLightTimeCalculator( 1, 0 );
    }

    LinkEnds getFirstLinkEnds( )
    {
        LinkEnds firstLinkEnds;
        firstLinkEnds[ transmitter ] = this->linkEnds_[ transmitter ];
        firstLinkEnds[ receiver ] = this->linkEnds_[ receiver ];
        return firstLinkEnds;
    }

    LinkEnds getSecondLinkEnds( )
    {
        LinkEnds secondLinkEnds;
        secondLinkEnds[ transmitter ] = this->linkEnds_[ transmitter2 ];
        secondLinkEnds[ receiver ] = this->linkEnds_[ receiver ];
        return secondLinkEnds;
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > >
    getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( 0, 0 ) } },
                 { std::make_pair( transmitter2, receiver ), { this->getSingleLegLightTimeCalculator( 1, 0 ) } } };
    }
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_POSITIONANGLEOBSERVATIONMODEL_H
