/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_HARMONIC_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_HARMONIC_GRAVITY_FIELD_VARIATIONS_H

#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/io/readHarmonicPotentialData.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace gravitation
{
enum class HarmonicGravityComponent { unspecified, oceanMass, directAtmosphericMass };
enum class HarmonicGravityLoading { unspecified, included, excluded };

//! Product-specific declarations, not inferred from potential amplitudes or a baseline field.
struct HarmonicGravityMetadata {
    HarmonicGravityComponent component = HarmonicGravityComponent::unspecified;
    double gravitationalParameter = TUDAT_NAN;
    double referenceRadius = TUDAT_NAN;
    bool fullyNormalized = false;
    std::string revision;
    HarmonicGravityLoading loading = HarmonicGravityLoading::unspecified;

    //! Archive the declared product conventions without reconstructing them from a file name.
    template< class Archive >
    void serialize( Archive& archive )
    {
        archive( component, gravitationalParameter, referenceRadius, fullyNormalized, revision, loading );
    }
};

//! Dimensionless normalized coefficients. Matrices use full degree/order indices.
//! One phase theta gives C=Ccos*cos(theta)+Csin*sin(theta), S=Scos*cos(theta)+Ssin*sin(theta).
struct HarmonicGravityConstituent {
    std::string doodson;
    std::array< int, 6 > multipliers{};
    Eigen::MatrixXd cCos, cSin, sCos, sSin;

    //! Archive the final phase identity and its four canonical coefficient matrices.
    template< class Archive >
    void serialize( Archive& archive )
    {
        archive( doodson, multipliers, cCos, cSin, sCos, sSin );
    }
};

//! Convert IERS 6.15 once; named profiles permit zeros only outside declared wave extents.
//! Generic sparse data rejects missing selected terms unless explicit zero-fill is requested.
//! Empty selectedDoodsonIds converts all in raw order; an explicit list preserves its supplied order.
std::vector< HarmonicGravityConstituent > convertHarmonicPotentialData( const input_output::HarmonicPotentialData& data,
                                                                        int maximumDegree,
                                                                        int maximumOrder,
                                                                        bool allowUnprofiledZeroFill = false,
                                                                        const std::vector< std::string >& selectedDoodsonIds = {} );

//! Exact ordered major->target map. Returned targets replace the major inventory, not append to it.
struct HarmonicGravityAdmittance {
    std::vector< std::string > majorDoodsonIds, targetDoodsonIds;
    Eigen::MatrixXd matrix;
};
std::vector< HarmonicGravityConstituent > applyHarmonicGravityAdmittance( const std::vector< HarmonicGravityConstituent >& majors,
                                                                          const HarmonicGravityAdmittance& admittance );

//! Common ocean/atmospheric potential evaluator. No extra phase, load or nodal multiplier.
class HarmonicGravityFieldVariations : public GravityFieldVariations
{
public:
    HarmonicGravityFieldVariations( const std::vector< HarmonicGravityConstituent >& constituents,
                                    const HarmonicGravityMetadata& metadata,
                                    const std::function< Eigen::Vector6d( double ) >& doodsonArguments,
                                    double outputGravitationalParameter,
                                    double outputReferenceRadius,
                                    int minimumDegree,
                                    int minimumOrder,
                                    int maximumDegree,
                                    int maximumOrder );
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( double time ) override;

private:
    std::vector< HarmonicGravityConstituent > constituents_;
    std::function< Eigen::Vector6d( double ) > doodsonArguments_;
};
}  // namespace gravitation
}  // namespace tudat
#endif
