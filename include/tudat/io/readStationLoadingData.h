#ifndef TUDAT_READ_STATION_LOADING_DATA_H
#define TUDAT_READ_STATION_LOADING_DATA_H

#include <istream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>
#include "tudat/basics/timeType.h"

namespace tudat
{
namespace input_output
{

//! Site of a station loading series, with its position and displacements.
struct StationLoadingSite {
    //! Cartesian position of the site (metres), and the distance (metres) from it within which a station may use its data.
    Eigen::Vector3d cartesianPosition = Eigen::Vector3d::Zero( );
    double validityRadiusMetres = 0.0;

    //! Geocentric radial, east and north displacements (metres), with the index of their epoch in StationLoadingData::epochs
    //! as key. Displacements may be missing at the start or end of the series, but not in between.
    std::map< std::size_t, Eigen::Vector3d > displacements;
};

//! Station loading series (for instance from the International Mass Loading Service): displacements of a set of sites at
//! a set of regularly spaced epochs.
struct StationLoadingData {
    //! Name of the file (or other source) of the series.
    std::string source;

    //! Epochs of the displacements (TAI), in increasing order.
    std::vector< Time > epochs;

    //! Sites, with their identifier in the file as key.
    std::map< std::string, StationLoadingSite > sites;
};

//! Read a station loading series in the EPHEDISP format (version 2005.06.30).
StationLoadingData readEphedispStationLoadingData( std::istream& stream, const std::string& source = "<stream>" );
StationLoadingData readEphedispStationLoadingData( const std::string& path );

}  // namespace input_output
}  // namespace tudat
#endif
