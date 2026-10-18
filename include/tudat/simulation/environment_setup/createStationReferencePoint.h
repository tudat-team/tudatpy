/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_CREATE_STATION_REFERENCE_POINT_H
#define TUDAT_CREATE_STATION_REFERENCE_POINT_H
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/ground_stations/stationReferencePoint.h"
namespace tudat
{
namespace simulation_setup
{
//! Create settings for a ground station at a reference point of a station (for instance its antenna reference point).
/*!
 *  The position of the station is that of the station coordinate history plus the reference point offsets. The station
 *  is named after the end point of the last offsets.
 *  \param history Station coordinate history (for instance of the station marker)
 *  \param referencePointOffsets Offset histories from the point of the history to the reference point, in order
 *  \param nominalPosition Nominal position of the ground station
 *  \return Ground station settings
 */
inline std::shared_ptr< GroundStationSettings > stationReferencePointSettings(
        std::shared_ptr< ground_stations::StationCoordinateHistory > history,
        std::vector< std::shared_ptr< ground_stations::StationReferencePointOffsetHistory > > referencePointOffsets,
        const Eigen::Vector3d& nominalPosition )
{
    if( !history || referencePointOffsets.empty( ) || !referencePointOffsets.back( ) )
    {
        throw std::runtime_error( "Reference-point station requires a station history and reference point offsets" );
    }
    const std::string name = referencePointOffsets.back( )->getOffsets( ).front( ).to.name;
    return std::make_shared< GroundStationSettings >(
            name,
            nominalPosition,
            coordinate_conversions::cartesian_position,
            std::vector< std::shared_ptr< GroundStationMotionSettings > >{
                    stationHistoryGroundStationMotionSettings( std::move( history ), std::move( referencePointOffsets ) ) } );
}
}  // namespace simulation_setup
}  // namespace tudat
#endif
