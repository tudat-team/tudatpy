/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_CREATE_STATION_REFERENCE_POINT_H
#define TUDAT_CREATE_STATION_REFERENCE_POINT_H
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/ground_stations/stationReferencePoint.h"
namespace tudat
{
namespace simulation_setup
{
//! Named observable point uses the selected absolute trajectory exactly once.
inline std::shared_ptr< GroundStationSettings > stationReferencePointSettings(
        std::shared_ptr< ground_stations::StationHistoryTrajectory > trajectory,
        const Eigen::Vector3d& nominalPosition,
        std::function< ground_stations::StationHistoryTime( const Time&, bool ) > sourceTimeFunction = nullptr )
{
    if( !trajectory || trajectory->getMetadata( ).referencePointKind == ground_stations::StationReferencePointKind::marker ||
        trajectory->getMetadata( ).referencePointName.empty( ) )
        throw std::runtime_error( "Reference-point station requires an explicit instrument point trajectory" );
    const std::string name = trajectory->getMetadata( ).referencePointName;
    return std::make_shared< GroundStationSettings >(
            name,
            nominalPosition,
            coordinate_conversions::cartesian_position,
            std::vector< std::shared_ptr< GroundStationMotionSettings > >{
                    stationHistoryGroundStationMotionSettings( std::move( trajectory ), std::move( sourceTimeFunction ) ) } );
}
}  // namespace simulation_setup
}  // namespace tudat
#endif
