#ifndef TUDAT_READ_ATMOSPHERIC_TIDE_DATA_H
#define TUDAT_READ_ATMOSPHERIC_TIDE_DATA_H

#include <istream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace tudat
{
namespace input_output
{

//! Coefficients (metres) of S1/S2 atmospheric tide loading at one point: rows radial, north and east displacement; columns
//! cos S1, sin S1, cos S2 and sin S2.
typedef Eigen::Matrix< double, 3, 4 > AtmosphericTideCoefficients;

//! Regular latitude/longitude grid of atmospheric tide loading coefficients (for instance the Ray and Ponte (2003) model).
struct AtmosphericTideGrid {
    //! Longitudes and latitudes (degrees) of the grid points.
    std::vector< double > longitudesDegrees;
    std::vector< double > latitudesDegrees;

    //! Coefficients at the grid points, longitude-major (index longitude * number of latitudes + latitude).
    std::vector< AtmosphericTideCoefficients > coefficientsMetres;
};

//! Read a grid in the Bernese 14-column format (longitude, latitude, 12 coefficients in millimetres; '$' header lines).
AtmosphericTideGrid readAtmosphericTideGrid( std::istream& stream, const std::string& source = "<stream>" );
AtmosphericTideGrid readAtmosphericTideGrid( const std::string& path );

//! Read station coefficients in the Bernese ATL format (station name, then three rows of four coefficients in millimetres;
//! '$' header lines), with the station name as key.
std::map< std::string, AtmosphericTideCoefficients > readAtmosphericTideStationData( std::istream& stream,
                                                                                     const std::string& source = "<stream>" );
std::map< std::string, AtmosphericTideCoefficients > readAtmosphericTideStationData( const std::string& path );

}  // namespace input_output
}  // namespace tudat
#endif
