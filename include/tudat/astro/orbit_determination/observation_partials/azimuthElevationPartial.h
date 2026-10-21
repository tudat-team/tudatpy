/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_AZIMUTHELEVATIONPARTIAL_H
#define TUDAT_AZIMUTHELEVATIONPARTIAL_H

#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"

namespace tudat
{

namespace observation_partials
{

Eigen::Matrix2d calculatePartialOfAzimuthElevationWrtAngularPosition( const Eigen::Vector2d& angularPosition,
                                                                      const Eigen::Matrix3d& rotationFromInertialToTopocentricFrame,
                                                                      const bool invertLineOfSight = false );

//! Derived class for scaling three-dimensional position partial to azimuth/elevation observable partial.
class AzimuthElevationScaling : public DirectPositionPartialScaling< 2 >
{
public:
    AzimuthElevationScaling( const std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator,
                             const observation_models::LinkEndType stationLinkEndType ):
        DirectPositionPartialScaling< 2 >( observation_models::azimuth_elevation_angle ),
        pointingAnglesCalculator_( pointingAnglesCalculator ), stationLinkEndType_( stationLinkEndType )
    {}

    ~AzimuthElevationScaling( ) {}

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    Eigen::Matrix< double, 2, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return referenceScalingFactor_ * getLinkEndSign( linkEndType );
    }

    Eigen::Matrix< double, 2, 3 > getFixedTimePositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        return scalingFactor_ * getLinkEndSign( linkEndType );
    }

    Eigen::Vector2d getLightTimePartialScalingFactor( )
    {
        return referenceLightTimeCorrectionScaling_;
    }

    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:
    double getLinkEndSign( const observation_models::LinkEndType linkEndType )
    {
        return ( linkEndType == observation_models::transmitter ) ? -1.0 : 1.0;
    }

    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator_;

    observation_models::LinkEndType stationLinkEndType_;

    Eigen::Matrix< double, 2, 3 > scalingFactor_;

    Eigen::Matrix< double, 2, 3 > referenceScalingFactor_;

    Eigen::Vector2d referenceLightTimeCorrectionScaling_;

    observation_models::LinkEndType currentLinkEndType_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_AZIMUTHELEVATIONPARTIAL_H
