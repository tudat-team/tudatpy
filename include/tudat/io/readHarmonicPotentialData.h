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
//! Identified provider inventories and Stokes units; generic IERS columns require their own unit declaration.
enum class HarmonicPotentialProfile { iersPotential, fes2004, fes2014bV1, costgEot11a };

//! Accept the provider's omitted leading zero (55.565) and full identifier (055.565).
//! Decode decimal Doodson identifiers and the verified TiME22 codes A=10, B=11, Z=-1.
std::array< int, 6 > parseDoodsonMultipliers( const std::string& identifier );

//! Published Cp, Sp, Cm, Sm in declared Stokes units; no phase or loading conversion here.
struct HarmonicPotentialCoefficient {
    //! Ordered Cp, Sp, Cm, Sm values before the declared Stokes unit scale is applied.
    std::array< double, 4 > values{};

    //! One-based original file line for source tracing.
    std::size_t sourceLine = 0;
};
//! One provider constituent, retaining sparse degree/order coverage and the original phase label.
struct HarmonicPotentialConstituent {
    //! Original identifier, canonical seven-character identifier and provider wave name.
    std::string originalDoodson, doodson, name;

    //! Multipliers of the six Doodson arguments in their published order.
    std::array< int, 6 > multipliers{};

    //! Values keyed by degree and order; missing source rows are not replaced with zero.
    std::map< std::pair< int, int >, HarmonicPotentialCoefficient > coefficients;

    //! Identified profile extent; generic/unidentified sparse input keeps these unspecified.
    int minimumDegree = -1, maximumDegree = -1, maximumOrder = -1;
};
//! Preserve ordered provider inventory and missing entries. GM/R are not guessed from filenames.
struct HarmonicPotentialData {
    //! Source identifier and original explanatory header lines.
    std::string source;

    std::vector< std::string > headerLines;

    //! Explicit provider profile controlling inventory, extent and unit checks.
    HarmonicPotentialProfile profile = HarmonicPotentialProfile::iersPotential;

    //! Multiplier converting published amplitudes to dimensionless Stokes coefficients.
    double stokesUnitScale = 0.0;

    //! Constituents in original source order, each retaining its own sparse coverage.
    std::vector< HarmonicPotentialConstituent > constituents;
};

//! Check phase identities, finite amplitudes and the selected profile's exact inventory and extents.
void validateHarmonicPotentialData( const HarmonicPotentialData& data );

//! Only normalized IERS potential columns are accepted, not water-height coefficients.
//! The identified COST-G EOT11a profile resolves its contradictory leftover 'in cm' header.
HarmonicPotentialData readHarmonicPotentialData( std::istream& stream,
                                                 const std::string& source = "<stream>",
                                                 HarmonicPotentialProfile profile = HarmonicPotentialProfile::iersPotential );
//! Open a potential file and retain its path as the source identifier.
HarmonicPotentialData readHarmonicPotentialData( const std::string& fileName,
                                                 HarmonicPotentialProfile profile = HarmonicPotentialProfile::iersPotential );
}  // namespace input_output
}  // namespace tudat
#endif
