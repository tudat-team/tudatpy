#include "Tudat/Astrodynamics/GroundStations/groundStation.h"


namespace tudat
{

namespace ground_stations
{

bool isTargetInView( const double time,
                     const Eigen::Vector3d targetRelativeState,
                     const boost::shared_ptr< PointingAnglesCalculator > pointingAngleCalculator,
                     const double minimumElevationAngle )
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
