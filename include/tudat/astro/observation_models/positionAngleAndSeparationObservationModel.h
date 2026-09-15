/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONANGLEANDSEPARATIONOBSERVATIONMODEL_H
#define TUDAT_POSITIONANGLEANDSEPARATIONOBSERVATIONMODEL_H

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <Eigen/Core>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/observationModel.h"

namespace tudat
{

namespace observation_models
{

template< typename ScalarType >
Eigen::Matrix< ScalarType, 2, 1 > calculatePositionAngleAndSeparation(
        const Eigen::Matrix< ScalarType, 3, 1 >& relativePositionFirstTransmitter,
        const Eigen::Matrix< ScalarType, 3, 1 >& relativePositionSecondTransmitter,
        const Eigen::Matrix< ScalarType, 3, 1 >& j2000NorthPoleDirection,
        const bool calculatePositionAngle = true )
{
    const ScalarType singularityTolerance = static_cast< ScalarType >( 100 ) * std::numeric_limits< ScalarType >::epsilon( );
    const ScalarType firstRange = relativePositionFirstTransmitter.norm( );
    const ScalarType secondRange = relativePositionSecondTransmitter.norm( );
    const ScalarType poleNorm = j2000NorthPoleDirection.norm( );
    if( firstRange <= singularityTolerance || secondRange <= singularityTolerance )
    {
        throw std::runtime_error( "Cannot calculate position angle and separation distance for a zero line-of-sight vector." );
    }
    if( poleNorm <= singularityTolerance )
    {
        throw std::runtime_error( "Cannot calculate position angle using a zero north-pole direction." );
    }

    const Eigen::Matrix< ScalarType, 3, 1 > firstLineOfSight = relativePositionFirstTransmitter / firstRange;
    const Eigen::Matrix< ScalarType, 3, 1 > secondLineOfSight = relativePositionSecondTransmitter / secondRange;
    const Eigen::Matrix< ScalarType, 3, 1 > northPoleDirection = j2000NorthPoleDirection / poleNorm;

    const ScalarType sineOfSeparation = firstLineOfSight.cross( secondLineOfSight ).norm( );
    if( sineOfSeparation <= singularityTolerance )
    {
        throw std::runtime_error(
                "Position angle is undefined and separation-distance partials are singular for coincident or antipodal lines of sight." );
    }

    ScalarType positionAngle = static_cast< ScalarType >( 0 );
    if( calculatePositionAngle )
    {
        Eigen::Matrix< ScalarType, 3, 1 > eastDirection = northPoleDirection.cross( firstLineOfSight );
        if( eastDirection.norm( ) <= singularityTolerance )
        {
            throw std::runtime_error( "Position angle is singular when the first line of sight is parallel to the ICRF/J2000 north pole." );
        }
        eastDirection.normalize( );
        const Eigen::Matrix< ScalarType, 3, 1 > northDirection = firstLineOfSight.cross( eastDirection );

        // Match the unnormalised right-ascension convention: return the principal value in [-pi, pi].
        positionAngle = std::atan2( secondLineOfSight.dot( eastDirection ), secondLineOfSight.dot( northDirection ) );
    }

    const ScalarType separationDistance = std::atan2( sineOfSeparation, firstLineOfSight.dot( secondLineOfSight ) );
    return ( Eigen::Matrix< ScalarType, 2, 1 >( ) << positionAngle, separationDistance ).finished( );
}

//! Class for simulating combined position-angle and angular-separation observables.
/*!
 *  Class for simulating combined position-angle and angular-separation observables, using light-time
 *  (with light-time corrections) to determine the states of the link ends (two transmitters and receiver).
 *  Returns a size-2 observable: [position_angle; angular_separation].
 *  Position angle is measured from ICRF/J2000 north through east and returned in [-pi, pi], matching the
 *  unnormalised right-ascension convention. The fixed reference pole is intended to become configurable in a future update.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class PositionAngleAndSeparationObservationModel : public ObservationModel< 2, ObservationScalarType, TimeType >
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
     *  \param j2000NorthPoleDirection Direction of the ICRF/J2000 north pole in the frame of the link-end states.
     *  \param calculatePositionAngle Whether to calculate position angle. This is false only for the separation-only wrapper.
     */
    PositionAngleAndSeparationObservationModel(
            const LinkEnds linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                    lightTimeCalculatorFirstTransmitter,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                    lightTimeCalculatorSecondTransmitter,
            const std::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = nullptr,
            const Eigen::Vector3d& j2000NorthPoleDirection = Eigen::Vector3d::UnitZ( ),
            const bool calculatePositionAngle = true ):
        ObservationModel< 2, ObservationScalarType, TimeType >(
                position_angle_and_separation,
                linkEnds,
                observationBiasCalculator,
                createFullLinkLightTimeCalculators( lightTimeCalculatorFirstTransmitter, lightTimeCalculatorSecondTransmitter ) ),
        j2000NorthPoleDirection_( j2000NorthPoleDirection ), calculatePositionAngle_( calculatePositionAngle )
    {}

    //! Destructor
    ~PositionAngleAndSeparationObservationModel( ) {}

    //! Function to compute an ideal position-angle and separation-distance observation at a given time, between two transmitters.
    /*!
     *  This function computes an ideal position-angle and separation-distance observation at a given time, between two transmitters.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated position-angle and separation-distance observable values as [PA; separation].
     */
    Eigen::Matrix< ObservationScalarType, 2, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettingsInput = nullptr ) override
    {
        // Check link end associated with input time and compute observable.
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error(
                    "Error when calculating position angle and separation distance observation, link end associated with time is not "
                    "receiver." );
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

        const Eigen::Matrix< ObservationScalarType, 2, 1 > positionAngleAndSeparation =
                calculatePositionAngleAndSeparation( relativeStateTransmitter1,
                                                     relativeStateTransmitter2,
                                                     j2000NorthPoleDirection_.template cast< ObservationScalarType >( ),
                                                     calculatePositionAngle_ );

        // Set link end times and states.
        linkEndTimes.clear( );
        linkEndStates.clear( );

        linkEndStates.push_back( firstLinkEndStates.at( 0 ) );
        linkEndStates.push_back( secondLinkEndStates.at( 0 ) );
        linkEndStates.push_back( firstLinkEndStates.at( 1 ) );

        linkEndTimes.push_back( firstLinkEndTimes.at( 0 ) );
        linkEndTimes.push_back( secondLinkEndTimes.at( 0 ) );
        linkEndTimes.push_back( firstLinkEndTimes.at( 1 ) );

        return positionAngleAndSeparation;
    }

    //! Function to get the object to calculate light time between first transmitter and receiver.
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return this->getSingleLegLightTimeCalculator( 0, 0 );
    }

    //! Function to get the object to calculate light time between second transmitter and receiver.
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

    Eigen::Vector3d getJ2000NorthPoleDirection( ) const
    {
        return j2000NorthPoleDirection_;
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > >
    getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( 0, 0 ) } },
                 { std::make_pair( transmitter2, receiver ), { this->getSingleLegLightTimeCalculator( 1, 0 ) } } };
    }

private:
    Eigen::Vector3d j2000NorthPoleDirection_;
    bool calculatePositionAngle_;
};

//! Class for simulating position angle observables, derived from the combined position-angle and separation model.
/*!
 *  Class for simulating position angle observables, using the PositionAngleAndSeparationObservationModel
 *  internally and extracting the position angle component (first element).
 *  The position angle is measured from ICRF/J2000 north through east at the first transmitter's line of sight.
 *  Its principal value is returned in [-pi, pi], matching the unnormalised right-ascension convention.
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
     *  \param j2000NorthPoleDirection Direction of the ICRF/J2000 north pole in the frame of the link-end states.
     */
    PositionAngleObservationModel( const LinkEnds linkEnds,
                                   const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                           lightTimeCalculatorFirstTransmitter,
                                   const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                           lightTimeCalculatorSecondTransmitter,
                                   const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
                                   const Eigen::Vector3d& j2000NorthPoleDirection = Eigen::Vector3d::UnitZ( ) ):
        ObservationModel< 1, ObservationScalarType, TimeType >(
                position_angle,
                linkEnds,
                observationBiasCalculator,
                createFullLinkLightTimeCalculators( lightTimeCalculatorFirstTransmitter, lightTimeCalculatorSecondTransmitter ) )
    {
        // Create the internal combined model with no bias (bias is handled at this level).
        positionAngleAndSeparationModel_ =
                std::make_shared< PositionAngleAndSeparationObservationModel< ObservationScalarType, TimeType > >(
                        linkEnds,
                        lightTimeCalculatorFirstTransmitter,
                        lightTimeCalculatorSecondTransmitter,
                        nullptr,
                        j2000NorthPoleDirection );
    }

    //! Destructor
    ~PositionAngleObservationModel( ) {}

    //! Function to compute ideal position angle observation at given time, between two transmitters.
    /*!
     *  This function computes the ideal position angle observation by delegating to the internal
     *  PositionAngleAndSeparationObservationModel and extracting the first element (position angle).
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid
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
        // Delegate computation to the internal combined model.
        Eigen::Matrix< ObservationScalarType, 2, 1 > positionAngleAndSeparationObservation =
                positionAngleAndSeparationModel_->computeIdealObservationsWithLinkEndData(
                        time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySettingsInput );

        // Extract position angle (first component)
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << positionAngleAndSeparationObservation( 0 ) ).finished( );
    }

    //! Function to get the object to calculate light time between first transmitter and receiver.
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return this->getSingleLegLightTimeCalculator( 0, 0 );
    }

    //! Function to get the object to calculate light time between second transmitter and receiver.
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

    Eigen::Vector3d getJ2000NorthPoleDirection( ) const
    {
        return positionAngleAndSeparationModel_->getJ2000NorthPoleDirection( );
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > >
    getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( 0, 0 ) } },
                 { std::make_pair( transmitter2, receiver ), { this->getSingleLegLightTimeCalculator( 1, 0 ) } } };
    }

private:
    //! Internal combined model that performs the full computation.
    std::shared_ptr< PositionAngleAndSeparationObservationModel< ObservationScalarType, TimeType > > positionAngleAndSeparationModel_;
};

//! Class for simulating angular separation distance observables, derived from the combined position-angle and separation model.
/*!
 *  Class for simulating angular separation distance observables, using the PositionAngleAndSeparationObservationModel
 *  internally and extracting the angular separation component (second element).
 *  The angular separation distance is the great-circle angle between the two transmitters as seen from the receiver.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class SeparationObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
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
    SeparationObservationModel( const LinkEnds linkEnds,
                                const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                        lightTimeCalculatorFirstTransmitter,
                                const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                                        lightTimeCalculatorSecondTransmitter,
                                const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >(
                separation_distance,
                linkEnds,
                observationBiasCalculator,
                createFullLinkLightTimeCalculators( lightTimeCalculatorFirstTransmitter, lightTimeCalculatorSecondTransmitter ) )
    {
        // Create the internal combined model with no bias (bias is handled at this level).
        positionAngleAndSeparationModel_ =
                std::make_shared< PositionAngleAndSeparationObservationModel< ObservationScalarType, TimeType > >(
                        linkEnds,
                        lightTimeCalculatorFirstTransmitter,
                        lightTimeCalculatorSecondTransmitter,
                        nullptr,
                        Eigen::Vector3d::UnitZ( ),
                        false );
    }

    //! Destructor
    ~SeparationObservationModel( ) {}

    //! Function to compute ideal angular separation distance observation at given time, between two transmitters.
    /*!
     *  This function computes the ideal angular separation distance observation by delegating to the internal
     *  PositionAngleAndSeparationObservationModel and extracting the second element (angular separation).
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid
     *  \param linkEndTimes List of times at each link end during observation (returned by reference).
     *  \param linkEndStates List of states at each link end during observation (returned by reference).
     *  \return Calculated angular separation distance observable value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettingsInput = nullptr ) override
    {
        // Delegate computation to the internal combined model.
        Eigen::Matrix< ObservationScalarType, 2, 1 > positionAngleAndSeparationObservation =
                positionAngleAndSeparationModel_->computeIdealObservationsWithLinkEndData(
                        time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySettingsInput );

        // Extract angular separation (second component)
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << positionAngleAndSeparationObservation( 1 ) ).finished( );
    }

    //! Function to get the object to calculate light time between first transmitter and receiver.
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getLightTimeCalculatorFirstTransmitter( )
    {
        return this->getSingleLegLightTimeCalculator( 0, 0 );
    }

    //! Function to get the object to calculate light time between second transmitter and receiver.
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

private:
    //! Internal combined model that performs the full computation.
    std::shared_ptr< PositionAngleAndSeparationObservationModel< ObservationScalarType, TimeType > > positionAngleAndSeparationModel_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_POSITIONANGLEANDSEPARATIONOBSERVATIONMODEL_H
