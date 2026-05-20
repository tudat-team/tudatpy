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
class AngularPositionObservationModel : public ObservationModel< 2, ObservationScalarType, TimeType >
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
            const std::shared_ptr< ObservationBias< 2 > > observationBiasCalculator = nullptr,
            const bool normalizeRightAscension = false ):
        ObservationModel< 2, ObservationScalarType, TimeType >(
                angular_position,
                linkEnds,
                observationBiasCalculator,
                std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >{
                        std::make_shared< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                                std::vector< std::shared_ptr<
                                        observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{ lightTimeCalculator },
                                std::make_shared< LightTimeConvergenceCriteria >( ),
                                false ) } ),
        normalizeRightAscension_( normalizeRightAscension )
    {
    }

    //! Destructor
    ~AngularPositionObservationModel( ) {}

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
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetingsInput = nullptr )
    {
        std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetings;
        std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator =
                getLightTimeCalculator( );
        this->setFrequencyProperties( time, linkEndAssociatedWithTime, lightTimeCalculator, ancillarySetingsInput, ancillarySetings );

        if( ancillarySetings != nullptr )
        {
            throw std::runtime_error( "Error, calling angular position observable with ancillary settings, but none are supported." );
        }

        this->getFullLinkLightTimeCalculatorFromBase( )->calculateLightTimeWithLinkEndsStates(
                time, linkEndAssociatedWithTime, linkEndTimes, linkEndStates, ancillarySetings );

        Eigen::Matrix< ObservationScalarType, 3, 1 > relativePosition =
                linkEndStates.at( 0 ).template cast< ObservationScalarType >( ).segment( 0, 3 ) -
                linkEndStates.at( 1 ).template cast< ObservationScalarType >( ).segment( 0, 3 );

        // Return observable
        double rightAscension = 2.0 *
                std::atan( relativePosition[ 1 ] /
                           ( std::sqrt( relativePosition[ 0 ] * relativePosition[ 0 ] + relativePosition[ 1 ] * relativePosition[ 1 ] ) +
                             relativePosition[ 0 ] ) );
        double declination = mathematical_constants::PI / 2.0 - std::acos( relativePosition[ 2 ] / relativePosition.norm( ) );
        //        return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << sphericalRelativeCoordinates.z( ),
        //                 mathematical_constants::PI / 2.0 - sphericalRelativeCoordinates.y( ) ).finished( );
        if( !normalizeRightAscension_ )
        {
            return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << rightAscension, declination ).finished( );
        }
        else
        {
            return ( Eigen::Matrix< ObservationScalarType, 2, 1 >( ) << rightAscension * std::cos( declination ), declination ).finished( );

        }
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

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > > getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( ) } } };
    }

    bool getNormalizeRightAscension( )
    {
        return normalizeRightAscension_;
    }


private:
    bool normalizeRightAscension_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_ANGULARPOSITIONOBSERVATIONMODEL_H
