/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SEPARATIONOBSERVATIONMODEL_H
#define TUDAT_SEPARATIONOBSERVATIONMODEL_H

#include <map>
#include <Eigen/Core>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/positionAngleAndSeparationObservationModel.h"

namespace tudat
{

namespace observation_models
{

inline double getSeparationScalingFactor( const observation_models::LinkEndType referenceLinkEnd,
                                          const std::vector< Eigen::Vector6d >& linkEndStates,
                                          const std::vector< double >& linkEndTimes,
                                          const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings,
                                          const bool isFirstPartial )
{
    return 1.0;
}

//! Class for simulating angular separation distance observables, derived from the combined PS model.
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
        // Create internal PS model with no bias (bias is handled at this level)
        psModel_ = std::make_shared< PositionAngleAndSeparationObservationModel< ObservationScalarType, TimeType > >(
                linkEnds, lightTimeCalculatorFirstTransmitter, lightTimeCalculatorSecondTransmitter, nullptr );
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
        // Delegate computation to the internal PS model
        Eigen::Matrix< ObservationScalarType, 2, 1 > psObs = psModel_->computeIdealObservationsWithLinkEndData(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySettingsInput );

        // Extract angular separation (second component)
        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << psObs( 1 ) ).finished( );
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
    //! Internal PS model that performs the full computation
    std::shared_ptr< PositionAngleAndSeparationObservationModel< ObservationScalarType, TimeType > > psModel_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_SEPARATIONOBSERVATIONMODEL_H
