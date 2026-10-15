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
//! Independent physical mass contributions used to detect overlapping constituent inventories.
enum class HarmonicGravityComponent { unspecified, oceanMass, directAtmosphericMass };

//! Declared source loading convention; the evaluator never applies an additional loading response.
enum class HarmonicGravityLoading { unspecified, included, excluded };

//! Product-specific declarations, not inferred from potential amplitudes or a baseline field.
struct HarmonicGravityMetadata {
    //! Physical ocean or direct atmospheric mass represented by the coefficients.
    HarmonicGravityComponent component = HarmonicGravityComponent::unspecified;

    //! Source gravitational parameter [m^3/s^2] and reference radius [m].
    double gravitationalParameter = TUDAT_NAN;

    double referenceRadius = TUDAT_NAN;

    //! Whether the supplied Stokes coefficients use full normalization.
    bool fullyNormalized = false;

    //! Product identity and revision, including any explicitly declared reference convention.
    std::string revision;

    //! Whether the source already contains the load-induced potential.
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
    //! Canonical identifier and multipliers of the six Doodson phase arguments.
    std::string doodson;

    std::array< int, 6 > multipliers{};

    //! Full degree/order matrices in Ccos, Csin, Scos, Ssin order; all values are dimensionless.
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
    //! Exact ordered input and output constituent identities.
    std::vector< std::string > majorDoodsonIds, targetDoodsonIds;

    //! Dimensionless matrix with target rows and major columns.
    Eigen::MatrixXd matrix;
};

//! Apply one explicit ordered map to all four quadratures, replacing the input inventory.
std::vector< HarmonicGravityConstituent > applyHarmonicGravityAdmittance( const std::vector< HarmonicGravityConstituent >& majors,
                                                                          const HarmonicGravityAdmittance& admittance );

//! Common ocean/atmospheric potential evaluator. No extra phase, load or nodal multiplier.
class HarmonicGravityFieldVariations : public GravityFieldVariations
{
public:
    //! Validate and retain the selected degree/order block, scaling source GM/radius to the consuming field once.
    HarmonicGravityFieldVariations( const std::vector< HarmonicGravityConstituent >& constituents,
                                    const HarmonicGravityMetadata& metadata,
                                    const std::function< Eigen::Vector6d( double ) >& doodsonArguments,
                                    double outputGravitationalParameter,
                                    double outputReferenceRadius,
                                    int minimumDegree,
                                    int minimumOrder,
                                    int maximumDegree,
                                    int maximumOrder );

    //! Evaluate the correction block at TDB seconds since J2000 using the supplied Doodson arguments [rad].
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( double time ) override;

private:
    //! Selected coefficient blocks after conversion to the consuming field's GM and radius.
    std::vector< HarmonicGravityConstituent > constituents_;

    //! Six phase arguments [rad] at each input TDB epoch [s since J2000].
    std::function< Eigen::Vector6d( double ) > doodsonArguments_;
};
}  // namespace gravitation
}  // namespace tudat
#endif
