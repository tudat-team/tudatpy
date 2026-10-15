/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_READ_OCEAN_POLE_GRAVITY_DATA_H
#define TUDAT_READ_OCEAN_POLE_GRAVITY_DATA_H

#include <complex>
#include <cstddef>
#include <istream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace tudat
{
namespace input_output
{
//! Dimensionless normalized integrals in IERS 6.23, not station-displacement coefficients.
struct OceanPoleGravityCoefficient {
    std::complex< double > a, b;
    std::size_t sourceLine = 0;

    //! Archive the published complex values and original source row.
    template< class Archive >
    void serialize( Archive& archive )
    {
        archive( a, b, sourceLine );
    }
};

//! Sparse published values: a missing entry remains distinct from a published zero.
struct OceanPoleGravityData {
    std::string source;
    std::vector< std::string > headerLines;
    std::map< std::pair< int, int >, OceanPoleGravityCoefficient > coefficients;
    int maximumDegree = 0;

    //! Archive sparse coverage and provenance without filling missing source terms.
    template< class Archive >
    void serialize( Archive& archive )
    {
        archive( source, headerLines, coefficients, maximumDegree );
    }
};

void validateOceanPoleGravityData( const OceanPoleGravityData& data );

//! Six columns: n m A(real) B(real) A(imaginary) B(imaginary). No decompression.
OceanPoleGravityData readOceanPoleGravityData( std::istream& stream, const std::string& source = "<stream>" );
OceanPoleGravityData readOceanPoleGravityData( const std::string& fileName );

}  // namespace input_output
}  // namespace tudat
#endif
