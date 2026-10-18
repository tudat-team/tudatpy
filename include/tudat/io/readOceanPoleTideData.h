#ifndef TUDAT_READ_OCEAN_POLE_TIDE_DATA_H
#define TUDAT_READ_OCEAN_POLE_TIDE_DATA_H

#include <cstddef>
#include <istream>
#include <string>
#include <vector>
#include <Eigen/Core>

namespace tudat
{
namespace input_output
{

//! Regular grid of IERS ocean pole tide loading coefficients, as in the file.
struct OceanPoleTideGrid {
    //! Number of grid points in longitude and latitude.
    std::size_t longitudeCount = 0;
    std::size_t latitudeCount = 0;

    //! Longitude and latitude (degrees) of the first grid point, and the grid spacing (degrees).
    double firstLongitudeDegrees = 0.0;
    double firstLatitudeDegrees = 0.0;
    double longitudeStepDegrees = 0.0;
    double latitudeStepDegrees = 0.0;

    //! Dimensionless coefficients per grid point, longitude varying fastest: rows radial, north, east; columns real and
    //! imaginary part.
    std::vector< Eigen::Matrix< double, 3, 2 > > coefficients;
};

//! Read the IERS eight-column plain-text format (header lines with the grid geometry, then longitude, latitude and the six
//! coefficients per grid point). No decompression is performed.
OceanPoleTideGrid readOceanPoleTideGrid( std::istream& stream, const std::string& source = "<stream>" );
OceanPoleTideGrid readOceanPoleTideGrid( const std::string& filePath );

}  // namespace input_output
}  // namespace tudat
#endif
