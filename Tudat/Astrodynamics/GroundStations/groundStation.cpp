/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/GroundStations/groundStation.h"


namespace tudat
{

namespace ground_stations
{

//! Function to check whether a target is visible from a ground station, based on minimum allowed elevation angle.
bool isTargetInView(
        const double time, const Eigen::Vector3d targetRelativeState,
        const std::shared_ptr< PointingAnglesCalculator > pointingAngleCalculator, const double minimumElevationAngle )
{
    double elevationAngle = pointingAngleCalculator->calculateElevationAngle( targetRelativeState, time );

    bool isStationVisible;
    if( elevationAngle > minimumElevationAngle )
    {
        isStationVisible = true;
    }
    else
    {
        isStationVisible = false;
    }
    return isStationVisible;
}

}

}
