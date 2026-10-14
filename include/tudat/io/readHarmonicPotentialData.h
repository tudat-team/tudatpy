/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_READ_HARMONIC_POTENTIAL_DATA_H
#define TUDAT_READ_HARMONIC_POTENTIAL_DATA_H

#include <array>
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
enum class HarmonicPotentialProfile { iersPotential, fes2004, fes2014bV1, costgEot11a };

//! Accept the provider's omitted leading zero (55.565) and full identifier (055.565).
std::array< int, 6 > parseDoodsonMultipliers( const std::string& identifier );

//! Published Cp, Sp, Cm, Sm in declared Stokes units; no phase or loading conversion here.
struct HarmonicPotentialCoefficient {
    std::array< double, 4 > values{};
    std::size_t sourceLine = 0;
};
struct HarmonicPotentialConstituent {
    std::string originalDoodson, doodson, name;
    std::array< int, 6 > multipliers{};
    std::map< std::pair< int, int >, HarmonicPotentialCoefficient > coefficients;
    //! Identified profile extent; generic/unidentified sparse input keeps these unspecified.
    int minimumDegree = -1, maximumDegree = -1, maximumOrder = -1;
};
//! Preserve ordered provider inventory and missing entries. GM/R are not guessed from filenames.
struct HarmonicPotentialData {
    std::string source;
    std::vector< std::string > headerLines;
    HarmonicPotentialProfile profile = HarmonicPotentialProfile::iersPotential;
    double stokesUnitScale = 0.0;
    std::vector< HarmonicPotentialConstituent > constituents;
};
void validateHarmonicPotentialData( const HarmonicPotentialData& data );
//! Only normalized IERS potential columns are accepted, not water-height coefficients.
//! The identified COST-G EOT11a profile resolves its contradictory leftover 'in cm' header.
HarmonicPotentialData readHarmonicPotentialData( std::istream& stream,
                                                 const std::string& source = "<stream>",
                                                 HarmonicPotentialProfile profile = HarmonicPotentialProfile::iersPotential );
HarmonicPotentialData readHarmonicPotentialData( const std::string& fileName,
                                                 HarmonicPotentialProfile profile = HarmonicPotentialProfile::iersPotential );
}  // namespace input_output
}  // namespace tudat
#endif
