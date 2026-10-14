#ifndef TUDAT_READ_OCEAN_POLE_TIDE_DATA_H
#define TUDAT_READ_OCEAN_POLE_TIDE_DATA_H

#include <cstddef>
#include <istream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include "tudat/io/readOceanLoadingData.h"

namespace tudat
{
namespace input_output
{
//! Neutral, dimensionless IERS ocean pole-load coefficients; rows radial/north/east, columns real/imaginary.
using OceanPoleTideCoefficients = Eigen::Matrix< double, 3, 2 >;

//! Uniform provider grid. Nodes are stored longitude-fast, then latitude, with no interpolation in the input layer.
struct OceanPoleTideGrid {
    std::size_t longitudeCount = 0;
    std::size_t latitudeCount = 0;
    double firstLongitudeDegrees = 0.0;
    double firstLatitudeDegrees = 0.0;
    double longitudeStepDegrees = 0.0;
    double latitudeStepDegrees = 0.0;
    std::vector< OceanPoleTideCoefficients > coefficients;
    std::string source;
    std::vector< std::string > sourceHeaders;
    OceanLoadingReferenceOrigin referenceOrigin = OceanLoadingReferenceOrigin::unspecified;
};

void validateOceanPoleTideGrid( const OceanPoleTideGrid& grid );

//! Read the IERS eight-column plain-text format, retaining headers and origin declaration. No decompression is performed.
OceanPoleTideGrid readOceanPoleTideGrid( std::istream& stream, const std::string& source = "<stream>" );
OceanPoleTideGrid readOceanPoleTideGrid( const std::string& filePath );

}  // namespace input_output
}  // namespace tudat
#endif
