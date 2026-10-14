/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_READ_ICGEM_GRAVITY_FIELD_H
#define TUDAT_READ_ICGEM_GRAVITY_FIELD_H
#include <array>
#include <istream>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include "tudat/astro/basic_astro/tideSystem.h"

namespace tudat
{
namespace input_output
{
enum class IcgemNormalization { fullyNormalized, unnormalized };
enum class IcgemCoefficientErrors { none, calibrated, formal, calibratedAndFormal };
//! Published static values and optional sigmaC/sigmaS pairs, in the declared normalization.
struct IcgemGravityCoefficient {
    double cosine = 0.0;
    double sine = 0.0;
    std::optional< std::array< double, 2 > > calibratedErrors;
    std::optional< std::array< double, 2 > > formalErrors;
    std::size_t sourceLine = 0;
    std::string originalRecord;
};
//! Neutral sparse static-field records. Missing entries remain distinct from published zeros.
struct IcgemGravityFieldData {
    std::string source;
    std::vector< std::string > headerLines;
    std::string productType;
    std::string modelName;
    std::string format = "icgem1.0";
    double gravitationalParameter = 0.0;
    double referenceRadius = 0.0;
    int maximumDegree = 0;
    IcgemCoefficientErrors errors = IcgemCoefficientErrors::none;
    //! These defaults are defined by the 2023 format, not inferred from coefficients.
    IcgemNormalization normalization = IcgemNormalization::fullyNormalized;
    basic_astrodynamics::TideSystem tideSystem = basic_astrodynamics::TideSystem::unknown;
    std::map< std::pair< int, int >, IcgemGravityCoefficient > coefficients;
};
void validateIcgemGravityFieldData( const IcgemGravityFieldData& data );
//! Read static gfc records. Temporal records reject; unknown keywords and trailing text are comments.
IcgemGravityFieldData readIcgemGravityField( std::istream& stream, const std::string& source = "<stream>" );
IcgemGravityFieldData readIcgemGravityField( const std::string& fileName );
}  // namespace input_output
}  // namespace tudat
#endif
