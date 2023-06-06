/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/corrections/solarCoronaCorrection.h"

namespace tudat
{

namespace observation_models
{

double SolarCoronaCorrection::computeMinimumDistanceOfLineOfSight(
        Eigen::Vector6d transmitterState,
        Eigen::Vector6d receiverState,
        Eigen::Vector6d sunState )
{
    Eigen::Vector3d sunPosition = sunState.segment( 0, 3 );
    Eigen::Vector3d transmitterRelativePosition = transmitterState.segment( 0, 3 ) - sunPosition;
    Eigen::Vector3d receiverRelativePosition = receiverState.segment( 0, 3 ) - sunPosition;

    // Moyer (2000), eq. 10-58
    Eigen::Vector3d lineOfSight = ( receiverRelativePosition - transmitterRelativePosition ).normalized( );

    // Moyer (2000), eq. 10-59
    Eigen::Vector3d minimumDistancePoint =
            transmitterRelativePosition - transmitterRelativePosition.dot( lineOfSight ) * lineOfSight;

    // Check if point is valid. If not, set minimum distance to NAN
    // Moyer (2000), eqs. 10-67, 10-68
    if ( transmitterRelativePosition.dot( lineOfSight ) < 0 && receiverRelativePosition.dot( lineOfSight ) > 0 )
    {
        return minimumDistancePoint.norm( );
    }
    else
    {
        return TUDAT_NAN;
    }

}


} // namespace observation_models

} // namespace tudat